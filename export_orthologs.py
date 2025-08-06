import argparse
import itertools
import sys
import re
import psutil
from find_homologs import highest_bitscores, eprint
from filtered_distance import (
    SampleSimilarity,
    get_ideal_components,
)
from path_to_sample import path_to_sample
from collections import defaultdict
from make_subset import multi_glob
from contextlib import ExitStack
from simple_blast import TabularBlastnSearch
from build_graph import component_subgraphs
from strand_sat import sat_assign_strands

import networkx as nx
import numpy as np

from pathlib import Path

import Bio.SeqIO
import Bio.Align

from tqdm import tqdm

from joblib import Parallel, delayed

try:
    import resource
except ImportError:
    resource = None
#from IPython import embed

default_gene_re = re.compile("^.*g([0-9]+)_i([0-9]+)")

def named_reverse_complement(t):
    rc = t.reverse_complement()
    rc.id = f"-{t.id}"
    rc.name = t.name
    rc.description = t.description
    return rc

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--graph", type=Path)
    parser.add_argument(
        "-c",
        "--comparisons",
        type=Path,
        nargs="+",
    )
    parser.add_argument(
        "-A",
        "--analysis-root",
        type=Path
    )
    parser.add_argument(
        "-s",
        "--samples",
        type=int
    )
    parser.add_argument(
        "--gene-regex",
        "-r",
        type=re.compile,
        help="Python regex for parsing sequence IDs",
        default=default_gene_re
    )
    parser.add_argument(
        "--out-dir",
        "-O",
        type=Path,
        required=True
    )
    parser.add_argument(
        "--by",
        "-b",
        choices=["sample", "component"],
        default="sample"
    )
    parser.add_argument(
        "--remove-non-contributing",
        "-N",
        action="store_true",
        help="Remove ideal components that do not contribute to the distance"
    )
    parser.add_argument(
        "--debug",
        action="store_true"
    )
    parser.add_argument(
        "-o",
        "--concat-id-order",
        choices=["before", "after"],
        default="after",
    )
    parser.add_argument(
        "--no-fix-strand",
        action="store_true"
    )
    parser.add_argument(
        "-i",
        "--allow-inconsistent",
        action="store_true"
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1
    )
    # parser.add_argument(
    #     "--collapse-duplicates",
    #     "-R",
    #     action="store_true",
    #     help="Collapse exact duplicate transcripts"
    # )
    args = parser.parse_args()
    if not args.graph:
        if args.analysis_root:
            args.graph = args.analysis_root / "graph.pkl"
        else:
            parser.error("Must provide either --graph or --analysis-root.")
    if not args.comparisons:
        if args.analysis_root:
            args.comparisons = list(
                multi_glob(args.analysis_root / "od2", ["*.pkl", "*.h5"])
            )
        else:
            parser.error(
                "Must provide either --comparisons or --analysis-root."
            )
    return args

def renamed_seqs(
        rename,
        s_seqs,
        attr="id",
        unset={"id", "title", "description"}
):
    unset = unset - {attr}
    for i, t in enumerate(s_seqs):
        setattr(t[-1], attr, rename(*t[:-1], getattr(t[-1], attr)))
        for s in unset:
            setattr(t[-1], s, "")
        # q.description = "ideal_component_{}".format(
        #     order[(s, g)]
        # )
        yield t

def seq_tuples(sample, gene_regex):
    for seq in Bio.SeqIO.parse(sample, "fasta"):
        match_ = gene_regex.search(seq.id)
        gene = int(match_.group(1))
        isoform = int(match_.group(2))
        yield (sample, gene, isoform, seq)

def concat_names(rename, order="after", sep=":"):
    def inner(sample_id, gene_id, isoform_id, old):
        new_component = rename(sample_id, gene_id, isoform_id)
        if order == "before":
            t = (new_component, old)
        else:
            t = (old, new_component)
        return sep.join(t)
    return inner

def get_strand(aligner, a, b):
    return 2*int(
        aligner.score(a, b) > aligner.score(a, b.reverse_complement())
    ) - 1

def is_mismatch(g, valid_genes, e):
    return all(t[:2] in valid_genes for t in e) and \
        g.nodes[e[0]]["strand"] != \
            g.nodes[e[1]]["strand"] * g.edges[e]["weight"]

def blast_pairwise_get_strands(isoforms):
    if len(isoforms) <= 1:
        return
    name_to_isoform = {i[1].id: i[0] for i in isoforms}
    seqs = [i[1] for i in isoforms]
    with TabularBlastnSearch.from_sequences(
            seqs,
            seqs,
            evalue=1e-5,
            additional_columns=["sstrand"]
    ) as search:
        hits = search.hits.loc[search.hits["qseqid"] > search.hits["sseqid"]]
        #print(hits[["qseqid", "sseqid"]])
        for _, x in highest_bitscores(
                hits,
                groupby=["qseqid", "sseqid"]
        ).iterrows():
            yield (
                (name_to_isoform[x["qseqid"]], x["qseqid"]),
                (name_to_isoform[x["sseqid"]], x["sseqid"]),
                2 * (x["sstrand"] == "plus") - 1
            )

def parallel_get_strands(gene_to_isoforms, index, jobs=1):
    return Parallel(n_jobs=jobs)(
        delayed(lambda x: list(blast_pairwise_get_strands(x)))(
            [(i[0], index[i[1]]) for i in isoforms]
        )
        for (gene, isoforms) in tqdm(gene_to_isoforms.items())
        if len(isoforms) > 1
    )

def build_strand_graph(sim, component_sample_genes, gene_regex, jobs=1):
    strand_graph = nx.Graph()
    # Add edges for isoform-isoform strands.
    for sample in sim.samples:
        index = Bio.SeqIO.index(sample, "fasta")
        gene_to_isoforms = defaultdict(list)
        for s in index:
            parsed = gene_regex.search(s)
            gene, isoform = map(int, parsed.groups())
            if (sample, gene) in component_sample_genes:
                gene_to_isoforms[gene].append((isoform, s))
        strand_graph.add_nodes_from(
            [
                (sample, gene, isoform)
                for (gene, isoforms) in gene_to_isoforms.items()
                for (isoform, _) in isoforms
            ]
        )
        strand_graph.add_weighted_edges_from(
            (
                (
                    sample,
                    gene,
                    ia,
                ),
                (
                    sample,
                    gene,
                    ib
                ),
                #get_strand(aligner, index[sa], index[sb]),
                strand
            )
            for it in parallel_get_strands(
                    gene_to_isoforms,
                    index,
                    jobs
            )
            for (ia, sa) , (ib, sb), strand in it
        )
    # Add edges for gene-gene strands.
    plus_minus = {"plus": 1, "minus": -1}
    plus_minus_inv = {v: k for (k, v) in plus_minus.items()}
    restricted_comparison_dfs = list(sim.restricted_comparison_dfs())
    for _, df in restricted_comparison_dfs:
        if "sstrand" not in df.columns:
            df["sstrand"] = np.sign(
                df["send"] - df["sstart"]
            ).apply(plus_minus_inv.__getitem__)
    strand_graph.add_weighted_edges_from(
        tuple(
            tuple(row[x + col] for col in ["sample", "gene", "iso"])
            for x in ["q", "s"]
        ) + (plus_minus[row["sstrand"]],)
        for _, df in restricted_comparison_dfs
        for _, row in df.iterrows()
    )
    components = component_subgraphs(strand_graph)
    gene_to_components = defaultdict(set)
    for component in components:
        for sample, gene, _  in component.nodes:
            gene_to_components[(sample, gene)].add(component)
    component_graph = nx.Graph()
    for gene_components in gene_to_components.values():
        component_graph.add_nodes_from(gene_components)
        component_graph.add_edges_from(
            itertools.pairwise(gene_components)
        )
    component_components = list(component_subgraphs(component_graph))
    node_to_component_component = {}
    for component_component in component_components:
        for subgraph in component_component.nodes:
            for node in subgraph.nodes:
                node_to_component_component[node] = component_component
            n1 = next(iter(subgraph.nodes))
            subgraph.nodes[n1]["strand"] = 1
            for a, b in nx.dfs_edges(subgraph, n1):
                subgraph.nodes[b]["strand"] = \
                    subgraph.nodes[a]["strand"] \
                    * subgraph.edges[(a, b)]["weight"]
    return strand_graph, node_to_component_component

def get_sample_gene_to_component(ideal):
    return {
        v: i for (i, c) in enumerate(ideal) for v in c.nodes
    }

class OrthologExporter:
    def __init__(
            self,
            sim,
            gene_regex=default_gene_re,
            non_contributing=True,
            collapse=False,
            consistent_strands=True,
            allow_inconsistent=True,
            jobs=1,
            debug=False,
    ):
        self.samples = sim.samples
        self.ideal = list(get_ideal_components(sim.graph, sim.sample_count))
        self.sample_gene_to_component = get_sample_gene_to_component(self.ideal)
        #print(self.sample_gene_to_component)
        self.gene_regex = gene_regex
        self.ideal_ids = set(range(len(self.ideal)))
        if not non_contributing:
            print("Filtering non-contributing.")
            total_distances = [0]*len(self.ideal)
            for samples, df in sim.restricted_comparison_dfs():
                for ix, row in df.iterrows():
                    dist = row["length"] - row["gaps"] - row["nident"]
                    assert dist >= 0
                    total_distances[
                        self.sample_gene_to_component[
                            row["qsample"],
                            row["qgene"]
                        ]
                    ] += dist
            if debug:
                for (i, d) in enumerate(total_distances):
                    if not d:
                        print(f"Excluding non-contributing component {i}.")
            self.ideal_ids = set(
                i for (i, d) in enumerate(total_distances) if d > 0
            )
            self.sample_gene_to_component = {
                k: v for (k, v) in self.sample_gene_to_component.items()
                if v in self.ideal_ids
            }
        if consistent_strands:
            print("Determining relative strand orientations.")
            # parse = functools.partial(parse_seq_id, self.gene_regex)
            # aligner = Bio.Align.PairwiseAligner()
            # aligner.substitution_matrix = Bio.Align.substitution_matrices.load(
            #     "NUC.4.4"
            # )
            # aligner.open_gap_score = -10
            # aligner.extend_gap_score = -0.5
            self.strand_graph, self.node_to_component_component = \
                build_strand_graph(
                    sim,
                    self.sample_gene_to_component,
                    self.gene_regex,
                    jobs=jobs
                )
            valid_genes = {tuple(x) for x in sim.valid.itertuples(index=False)}
            mismatches = {
                    e
                    for e in self.strand_graph.edges
                    if is_mismatch(self.strand_graph, valid_genes, e)
            }
            mismatch_component_components= {
                self.node_to_component_component[m[0]]
                for m in mismatches
            }
            # mismatch_components = list(
            #     map(
            #         strand_graph.subgraph,
            #         mismatch_component_nodes
            #     )
            # )
            for m in mismatches:
                self.strand_graph.edges[m]["mismatch"] = True
            if mismatches:
                eprint(
                    "Found {} strand mismatches in {} component groups.".format(
                        len(mismatches),
                        len(mismatch_component_components)
                    )
                )
                if not allow_inconsistent:
                    sys.exit(1)
                else:
                    eprint("Attempting fix.")
                for i, comp_comp in enumerate(mismatch_component_components):
                    subgraph = self.strand_graph.subgraph(
                        n for comp in comp_comp for n in comp.nodes
                    )
                    # mm = sum(
                    #     1
                    #     for e in subgraph.edges
                    #     if is_mismatch(strand_graph, valid_genes, e)
                    # )
                    cost = sat_assign_strands(subgraph)
                    mm2 = sum(
                        1
                        for e in subgraph.edges
                        if is_mismatch(self.strand_graph, valid_genes, e)
                    )
                    if cost != mm2:
                        print("Bad cost!")
                        from IPython import embed; embed()
                    # print(
                    #     "Component {} reassignment: {} -> {}".format(
                    #         i,
                    #         mm,
                    #         mm2
                    #     )
                    # )
                mm2 = sum(
                    1
                    for e in self.strand_graph.edges
                    if is_mismatch(self.strand_graph, valid_genes, e)
                )
                eprint(
                    "Reduced mismatches: {} -> {}".format(len(mismatches), mm2)
                )
        #self.non_contributing = non_contributing
        self.collapse = collapse
        #self.valid_tuples = set(map(tuple, sim.valid.itertuples(index=False)))
        #self.samples = sim.valid["sample"].drop_duplicates()

    def _name_ideal(self, sample, gene, isoform):
        return "ideal_component_{}".format(
            self.sample_gene_to_component[(sample, gene)]
        )

    def _orient(self, t):
        try:
            if self.strand_graph.nodes[t[:-1]]["strand"] == -1:
                return t[:-1] + (
                    named_reverse_complement(
                        t[-1]
                    ),
                )                    
        except AttributeError:
            pass
        return t

    def by_sample(self, out_dir, rename=None, order="after"):
        if rename is None:
            rename = concat_names(self._name_ideal, order=order)
        for sample in self.samples:
            out_fn = out_dir / "{}_orthologs.fasta".format(
                path_to_sample(sample)
            )
            Bio.SeqIO.write(
                (
                    t[-1] for t in
                    renamed_seqs(
                        rename,
                        sorted(
                            (
                                self._orient(t)
                                for t in seq_tuples(sample, self.gene_regex)
                                if self.sample_gene_to_component.get(t[:-1])
                                in self.ideal_ids
                            ),
                            key=lambda x: self.sample_gene_to_component[x[:-1]]
                        ),
                    )
                ),
                out_fn,
                "fasta"
            )

    def by_component(
            self,
            out_dir,
            rename=None,
            order="after",
            set_rlimit=True
    ):
        if rename is None:
            rename = concat_names(
                lambda a, b, c: path_to_sample(a),
                order=order
            )
        component_paths = {
            i: out_dir / f"ideal_component_{i}.fasta"
            for i in self.ideal_ids
        }
        if set_rlimit:
            try:
                file_max = resource.getrlimit(resource.RLIMIT_OFILE)[1]
                resource.setrlimit(
                    resource.RLIMIT_OFILE,
                    (
                        min(
                            (len(component_paths) + \
                             len(psutil.Process().open_files())) * 2,
                            file_max
                        ),
                        file_max
                    )
                )
                print("rlimit is", resource.getrlimit(resource.RLIMIT_OFILE))
            except AttributeError:
                from IPython import embed; embed()
            
        #print(component_paths)
        with ExitStack() as stack:
            component_files = {
                i: stack.enter_context(open(f, "w"))
                for (i, f) in component_paths.items()
            }
            for sample in self.samples:
                # for _, gene, isoform, seq in seq_tuples(
                #         sample,
                #         self.gene_regex
                # ):
                #     pass
                print(sample)
                #from IPython import embed; embed()
                for _, gene, isoform, seq in renamed_seqs(
                        rename,
                        (
                            t for t in seq_tuples(sample, self.gene_regex)
                            if t[:-2] in self.sample_gene_to_component
                        )
                ):
                    Bio.SeqIO.write(
                        self._orient((sample, gene, isoform, seq))[-1],
                        component_files[
                            self.sample_gene_to_component[
                                (
                                    sample,
                                    int(
                                        self.gene_regex.search(seq.id).group(1)
                                    )
                                )
                            ]
                        ],
                        "fasta"
                    )
        return component_paths
        # if self.collapse:
        #     for f in component_paths.values():
        #         duplicates = defaultdict(list)
        #         for seq in Bio.SeqIO.parse(f, "r") as seq:
        #             duplicates[seq.seq].append(seq.description)

def main():
    args = handle_arguments()
    args.out_dir.mkdir(exist_ok=True)
    sim = SampleSimilarity.from_filenames(
        args.graph,
        tqdm(args.comparisons),
        sample_count=args.samples,
        store_dfs=True
    )
    exporter = OrthologExporter(
        sim,
        args.gene_regex,
        not args.remove_non_contributing,
        debug=args.debug,
        consistent_strands=not args.no_fix_strand,
        allow_inconsistent=args.allow_inconsistent,
        jobs=args.jobs,
    )
    getattr(exporter, "by_{}".format(args.by))(
        args.out_dir,
        order=args.concat_id_order
    )


if __name__ == "__main__":
    main()
