import argparse
import pickle
import functools
import re
from filtered_distance import SampleSimilarity, get_ideal_components
from path_to_sample import path_to_sample
from collections import Counter, defaultdict
from make_subset import multi_glob
from contextlib import ExitStack

from pathlib import Path

import Bio.SeqIO

from tqdm import tqdm
#from IPython import embed

default_gene_re = re.compile("^.*g([0-9]+)_i([0-9]+)")

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
        yield t[-1]

def seq_tuples(sample, gene_regex):
    for seq in Bio.SeqIO.parse(sample, "fasta"):
        yield (
            sample,
            int(
                gene_regex.search(seq.id).group(1)
            ),
            seq
        )

def concat_names(rename, order="after"):
    def inner(sample_id, gene_id, old):
        new_component = rename(sample_id, gene_id)
        if order == "before":
            t = (new_component, old)
        else:
            t = (old, new_component)            
        return "{} {}".format(*t)
    return inner

    

class OrthologExporter:
    def __init__(
            self,
            sim,
            gene_regex=default_gene_re,
            non_contributing=True,
            collapse=False,
            debug=False,
    ):
        self.samples = sim.samples
        self.ideal = list(get_ideal_components(sim.graph, sim.sample_count))
        self.sample_gene_to_component = {
            v: i for (i, c) in enumerate(self.ideal) for v in c.nodes
        }
        self.gene_regex = gene_regex
        self.ideal_ids = set(range(len(self.ideal)))
        if not non_contributing:
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
        #self.non_contributing = non_contributing
        self.collapse = collapse
        #self.valid_tuples = set(map(tuple, sim.valid.itertuples(index=False)))
        #self.samples = sim.valid["sample"].drop_duplicates()

    def _name_ideal(self, *t):
        return "ideal_component_{}".format(
            self.sample_gene_to_component[t]
        )


    def by_sample(self, out_dir, rename=None, order="after"):
        if rename is None:
            rename = concat_names(self._name_ideal, order=order)
        for sample in self.samples:
            out_fn = out_dir / "{}_orthologs.fasta".format(
                path_to_sample(sample)
            )
            Bio.SeqIO.write(
                renamed_seqs(
                    rename,
                    sorted(
                        (
                            t for t in seq_tuples(sample, self.gene_regex)
                            if self.sample_gene_to_component.get(t[:-1])
                            in self.ideal_ids
                        ),
                        key=lambda x: self.sample_gene_to_component[x[:-1]]
                    ),
                ),
                out_fn,
                "fasta"
            )

    def by_component(self, out_dir, rename=None, order="after"):
        if rename is None:
            rename = concat_names(
                lambda a, b: path_to_sample(a),
                order=order
            )
        component_paths = {
            i: out_dir / f"ideal_component_{i}.fasta"
            for i in self.ideal_ids
        }
        #print(component_paths)
        with ExitStack() as stack:
            component_files = {
                i: stack.enter_context(open(f, "w"))
                for (i, f) in component_paths.items()
            }
            for sample in self.samples:
                for seq in renamed_seqs(
                        rename,
                        (
                            t for t in seq_tuples(sample, self.gene_regex)
                            if t[:-1] in self.sample_gene_to_component
                        )
                ):
                    Bio.SeqIO.write(
                        seq,
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
        debug=args.debug
    )
    getattr(exporter, "by_{}".format(args.by))(
        args.out_dir,
        order=args.concat_id_order
    )

    
if __name__ == "__main__":
    main()
