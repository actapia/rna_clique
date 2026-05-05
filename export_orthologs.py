import itertools
import re
try:
    import resource
except ImportError:
    resource = None

import psutil
import Bio
import Bio.SeqIO
import Bio.Align
import networkx as nx
import numpy as np

import config as config_module

from collections import defaultdict
from contextlib import ExitStack
from typing import Callable, Iterator, Optional, TypeVarTuple, Any
from collections.abc import Iterable, Mapping
from pathlib import Path

from simple_blast import TabularBlastnSearch
from tqdm import tqdm
from joblib import Parallel, delayed

from find_homologs import highest_bitscores
from filtered_distance import (
    SampleSimilarity,
    get_ideal_components,
)
from path_to_sample import path_to_sample
from graph import component_subgraphs
from strand_sat import sat_assign_strands
from transcripts import TranscriptID
from gene_matches_tables import get_table_files

default_gene_re = re.compile("^.*g([0-9]+)_i([0-9]+)")

def named_reverse_complement(t: Bio.SeqRecord) -> Bio.SeqRecord:
    """Get reverse complement sequence record with ID derived from original.

    This function returns the reverse complement of a Bio.SeqRecord object with
    the sequence ID changed to be the same as the original sequeunce, but with a
    "-" prepended.

    Parameters:
        t: The Bio.SeqRecord object for which to obtain the reverse complement.

    Returns:
        Reverse complement sequence record with ID derived from the original.
    """
    rc = t.reverse_complement()
    rc.id = f"-{t.id}"
    rc.name = t.name
    rc.description = t.description
    return rc

def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager(
        description="Export ortholog sequences from ideal components."
    )
    arg_config.expose_fields_with_default_aliases(
        "graph",
        "tables_dir",
        "jobs",
        "transcript_id_regex",
        required=True
    )
    arg_config.expose_config_field(
        "output_dir",
        aliases=["--analysis-root", "--rna-clique-output-dir", "-A"],
        help="RNA-clique analysis root (output_dir).",
    )
    arg_config.add_argument(
        "--export-output-dir",
        "-X",
        type=Path,
        required=True,
        help="Output directory in which to store exported orthologs.",
    )
    arg_config.add_argument(
        "--by",
        "-b",
        choices=["sample", "component"],
        default="sample",
        help="Attribute by which to organize orthologs in export.",
    )
    arg_config.add_argument(
        "--remove-non-contributing",
        "-N",
        action="store_true",
        help="Remove ideal components that contribute no differences.",
    )
    arg_config.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug behavior.",
    )
    arg_config.add_argument(
        "-o",
        "--concat-id-order",
        choices=["before", "after"],
        default="after",
        help="Where to place original sequence ID relative to group name.",
    )
    arg_config.add_argument(
        "--no-fix-strand",
        action="store_true",
        help="Do not attempt to put transcripts in consistent orientations.",
    )
    arg_config.add_argument(
        "-i",
        "--allow-inconsistent",
        action="store_true",
        help="Approximate transcript reorientation instead of failing.",
    )
    arg_config.add_argument(
        "--all",
        "-a",
        action="store_true",
        help="Create combined all_ideal.fasta file."
    )
    return arg_config

Ts = TypeVarTuple("Ts")

def renamed_seqs(
        rename: Callable[tuple[*Ts], str],
        s_seqs: Iterable[tuple[*Ts, Bio.SeqRecord]],
        attr: str = "id",
        unset: set[str] = {"id", "title", "description"}
) -> Iterator[tuple[*Ts, Bio.SeqRecord]]:
    """Yield sequence tuples renamed using the given function.

    The sequence tuples provided must contain information used by the provided
    function to generate the name (for example, the sample, gene, and sequence
    IDs) followed by the actual SeqRecord for the sequence. The SeqRecord itself
    is not passed to the provided function.

    The provided function should return a string, the new name for the sequence.

    Since SeqRecord has default values for certain attributes if they are unset,
    this function can also automatically set those values to the empty string,
    "".

    Parameters:
        rename:      Function to generate name from sequence tuple information.
        s_seqs:      Sequence tuples to be renamed.
        attr (str):  Attribute to set to "rename" a sequence.
        unset (set): Attributes to set to the empty string by default.
    """
    unset = unset - {attr}
    for i, t in enumerate(s_seqs):
        setattr(t[-1], attr, rename(*t[:-1], getattr(t[-1], attr)))
        for s in unset:
            setattr(t[-1], s, "")
        # q.description = "ideal_component_{}".format(
        #     order[(s, g)]
        # )
        yield t

def seq_tuples(
        sample: str | Path,
        parse_transcript_id: Callable[[str], TranscriptID]
) -> Iterator[tuple[str | Path, int, int, Bio.SeqRecord]]:
    """Iterate over tuples of sample path, gene, isoform, and sequence in file.

    Parameters:
        sample:              Path to FASTA file containing sample transcripts.
        parse_transcript_id: Function to parse transcript FASTA IDs.
    """
    for seq in Bio.SeqIO.parse(sample, "fasta"):
        _, gene, isoform = parse_transcript_id(seq.id)
        yield (sample, gene, isoform, seq)

def concat_names(
        rename: Callable[[str, int, int], str],
        order: str = "after",
        sep: str = ":"
) -> Callable[[str, int, int, str], str]:
    """Get a function that combines original name with new name from function.

    This function takes an existing function for generating a name for a
    sequence based on its sample name, gene ID, and isoform ID and returns a new
    function that generates a name that is the original name join with the new
    name produced by the original function.

    The first argument, rename, must be a function accepting three
    arguments. The first should be the sample ID (path or string representing
    the sample). The second and third should be the gene and isoform IDs of the
    sequence, respectively.

    The second argument, order, indicates the order in which the old name and
    new name should be joined. "before" indicates that the new name comes first;
    "after" indicates that the old name comes first.

    The third argument, sep, is the separator used to join the names.

    Parameters:
        rename:      Function mapping sample, transcript, and isoform ID to name
        order (str): Whether to put the new name "before" or "after".
        sep (str):   Separator for joining the old and new names.

    Returns:
        Function mapping sample, transcript, isoform IDs & old name to new name.
    """
    def inner(
            sample_id: str,
            gene_id: int,
            isoform_id: int,
            old: str
    ) -> str:
        new_component = rename(sample_id, gene_id, isoform_id)
        if order == "before":
            t = (new_component, old)
        else:
            t = (old, new_component)
        return sep.join(t)
    return inner

def get_strand(
        aligner: Bio.Align.PairwiseAligner,
        a: Bio.SeqRecord,
        b: Bio.SeqRecord
) -> int:
    """Use a PairwiseAligner to get the relative orientations of two sequences.

    Parameters:
        aligner: The PairwiseAligner to use for checking relative orientation.
        a:       The first sequence to align.
        b:       The second sequnece to align.

    Returns:
        1 if the two are in the same orientation, -1 otherwise.
    """
    return 2*int(
        aligner.score(a, b) > aligner.score(a, b.reverse_complement())
    ) - 1

def is_mismatch(
        g: nx.Graph,
        valid_genes: set[tuple[str, int]],
        e: tuple[str, int, int]
) -> bool:
    """Determine if an edge in an assigned graph is a "mismatch" edge.

    A mismatch edge is one whose endpoints are present in ideal components
    (i.e., are "valid") and for which the product of the relative orientation
    for the edge and one endpoint's assignment is not the same as the other
    other endpoint's assignment.

    This happens in exactly two ways. First, if two endpoints from ideal
    components are both assigned 1, but the edge corresponds to relative
    orientation -1, the edge is a mismatch. Second, if two endpoints from ideal
    components are assigned 1 and -1, but the edge is 1, then the edge is a
    mismatch.

    Parameters:
        g:                 Assigned graph containing the given edge.
        valid_genes (set): Set of genes in ideal components.
        e (tuple):         Edge checked for being a mismatch edge.

    Returns:
        Whether the provided edge e is a mismatch edge in the assigned graph g.
    """
    return all(t[:2] in valid_genes for t in e) and \
        g.nodes[e[0]]["strand"] != \
            g.nodes[e[1]]["strand"] * g.edges[e]["weight"]

def blast_pairwise_get_strands(
        isoforms: Iterable[tuple[int, Bio.SeqRecord]]
) -> Iterator[tuple[tuple[int, str], tuple[int, str], int]]:
    """Get the relative orientations of isotigs using BLAST.

    This function accepts an iterable of pairs of isoform IDs and their
    corresponding Bio.SeqRecords for isoforms of the same gene and yields
    3-tuples representing the relative orientations of pairs of input
    transcripts. In each tuple, the first and second elements are both pairs
    representing the transcripts for which the relative orientation is
    given. The first element of each of these tuples is the isoform ID of the
    isforom, and the second is the isotig FASTA sequence ID. The last element of
    the yielded tuple is an integer indicating the relative orientation of the
    isoforms represented by the first two elements of the tuple. This last
    element is +1 if the two isoforms are in the same orientation and is -1 if
    the two isoforms are in reverse complement orientation.

    Although this function usually obtains the correct relative orientations for
    all isoforms, it is ocassionally incorrect. Hence, some care should be taken
    when relying on this function.

    Parameters:
        isoforms: Itereable of isoform ID and corresponding SeqRecord pairs.
    """
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

def parallel_get_strands(
        gene_to_isoforms: Mapping[int, Iterable[tuple[int, str]]],
        index: Mapping[[str], Bio.SeqRecord],
        jobs: int = 1
) -> Iterable[list[tuple[tuple[int, str], tuple[int, str], int]]]:
    """Get pairwise relative orientations within sets of isoforms in parallel.

    This function uses the blast_pairwise_get_strands function to obtain
    pairwise relative orientations within multiple sets of isoforms.

    The function accepts a mapping from gene IDs to lists of pairs representing
    isoforms belonging to that gene. Each pair should consist of the isoform ID,
    followed by the FASTA sequence ID for that isoform.

    This function also requires an index, which should map FASTA sequence IDs of
    isoforms to SeqRecord objects. Such an index can be created using
    Bio.SeqIO.index.

    This function returns a rather complicated object---it is an Iterable of
    lists. Each list contains all elements yielded from
    blast_pairwise_get_strands for a single gene (isotig set); these are the
    relative orientations for the set of isoforms. Each element of such a list
    is a 3-tuple. The first two elements of each 3-tuple are both pairs
    representing the isoforms for which the relative orientation was found. The
    first element represents the isoform ID, and the second element is the FASTA
    sequence ID of the isoform. The third element of each 3-tuple is an integer
    indicating the relative orientation of the isoforms represented by the first
    two elements of the 3-tuple. This third element is +1 if the two isoforms
    are in the same orientation, and it is -1 if the two isoforms are in reverse
    complement orientation.

    Parameters:
        gene_to_isoforms (dict): Mapping from gene IDs to lists of isoforms.
        index:                   Mapping to retrieve SeqRecords from FASTA IDs.
        jobs (int):              Number of parallel jobs to use.

    Returns:
        Pairwise relative orientations for all provided isotigs sets.
    """
    return Parallel(n_jobs=jobs)(
        delayed(lambda x: list(blast_pairwise_get_strands(x)))(
            [(i[0], index[i[1]]) for i in isoforms]
        )
        for (gene, isoforms) in tqdm(gene_to_isoforms.items())
        if len(isoforms) > 1
    )

def build_strand_graph(
        sim: SampleSimilarity,
        component_sample_genes: dict[tuple[str, int], int],
        parse_transcript_id: Callable[[str], TranscriptID],
        jobs: int = 1
) -> tuple[nx.Graph, dict[tuple[str, int, int], nx.Graph]]:
    """Construct a graph representing relative transcript orientations.

    The strand graph consists of vertices representing specific isotigs; they
    are tuples of sample, gene, and isoform IDs. An edge exists between two
    isotigs when the relative orientation of the isotigs---either forward or
    reverse compelement---is known. When the two isotigs are in the same
    orientation, the edge has "weight" attribute +1. Otherwise, the edge has
    weight -1.

    This function also constructs a "meta-strand graph" in which each vertex is
    a connected component of the strand graph. There is a one-to-one
    correspondance between components of the meta-strand graph and ideal
    components in the gene matches graph.

    Detection of relative orientations is error prone; some edges may have
    inaccurate weights.

    Vertices in the strand graph can be assigned "strand" attributes indicating
    whether they need to be reoriented to ensure all isotigs in the component
    are in the same orientation. When the strand attribute is +1, no
    reorientation is necessary. When the strand attribute is -1, the isotig
    needs to be reoriented as its reverse complement. This function returns an
    unassigned graph; other functions can be used to assign the strand
    attributes based on the edge weights.

    This function also returns a dict mapping nodes of the strand graph to the
    components of the meta-strand graph to which they belong.

    The function requires the SampleSimilarity object for the analysis for which
    the strand graph should be built and a dictionary mapping sample and gene ID
    tuples to ideal component IDs.

    Parameters:
        sim:                           SampleSimilarity for the analysis
        component_sample_genes (dict): dict mapping (sample ID, gene ID) tuples
                                       to ideal component IDs.
        parse_transcript_id:           Function to parse transcript IDs from
                                       FASTA sequence IDs.
        jobs (int):                    Number of parallel jobs to use.

    Returns:
        Strand graph and dict from sample, gene, isoform IDs to meta-components.
    """
    strand_graph = nx.Graph()
    # Add edges for isoform-isoform strands.
    for sample in sim.samples:
        index = Bio.SeqIO.index(sample, "fasta")
        gene_to_isoforms = defaultdict(list)
        for s in index:
            _, gene, isoform = parse_transcript_id(s)
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
    return strand_graph, node_to_component_component

def dfs_assign_strands(strand_graph: nx.Graph):
    """Assign orientations to strand graph nodes using a depth-first traversal.
    
    This function assigns the "strand" attribute of strand graph nodes by first
    assigning 1 to an arbitrary node in the graph, and then beginning a
    depth-first traversal of the graph rooted at the initial arbitrary node. For
    each edge in the traversal, the successor node is assigned the product of
    the predecessor node's strand and the edge's weight, ensuring that the
    traversed edge is consistent.

    This function always produces a consistent assignment if one exists. That
    is, if there is an assignment of strands to nodes such that, for every edge
    in the graph, the product of the edge weight and the strand of one node is
    equal to the strand of the other node, then this function will find such an
    assignment.

    If there is no such assignment, this function can perform very poorly,
    producing many more edges for which the desired property is violated than
    necessary.

    Parameters:
        strand_graph: The strand graph for which to assign orientations.
    """
    components = list(component_subgraphs(strand_graph))
    for subgraph in components:
        n1 = next(iter(subgraph.nodes))
        subgraph.nodes[n1]["strand"] = 1
        for a, b in nx.dfs_edges(subgraph, n1):
            subgraph.nodes[b]["strand"] = \
                subgraph.nodes[a]["strand"] \
                * subgraph.edges[(a, b)]["weight"]


def get_sample_gene_to_component(
        ideal: list[nx.Graph]
) -> dict[tuple[str, int], int]:
    """Get a dict mapping sample-gene pairs to ideal component indices.

    The values of the dict returned are indices in the provided ideal component
    list. They are not the components themselves.

    Parameters:
        ideal (list): List of ideal components from a gene matches graph.

    Returns:
        A dict mapping each sample-gene pair to the index of its component.
    """
    return {
        v: i for (i, c) in enumerate(ideal) for v in c.nodes
    }

class InconsistentGraphError(ValueError):
    """Raised when a strand graph does not have a consistent assignment.

    A consistent assignment is one in which, for every edge in the graph, the
    product of the edge's weight (relative orientation) and one vertex's
    assignment is always the same as the other vertex's assignment.
    """
    pass

class OrthologExporter:
    """Class for exporting orthologs detected by RNA-clique.

    This class exports the sequences of transcripts belonging to genes in ideal
    components detected by RNA-clique. These include all isoforms of these
    genes---not merely those that actually appear in the ideal components.

    The class includes various options for controlling the way in which the
    transcripts are exported. First, the exporter can "slice" the orthologs two
    different ways by using either the by_sample or by_component method. The
    first groups sequences into FASTA files based on sample. The second groups
    seqsuences into FASTA files based on ideal component. In both cases, the
    ideal component and sample to which a sequence belongs can always be
    recovered, but certain ways of slicing the data will be preferred for
    certain applications.

    Second, the exporter can optionally exclude "non-contributing"
    components. Despite their name, these components do affect the
    distance. They are simply the components where there is no variation among
    the samples, and, thus, they do not contribute differences. For exports
    where the goal is to inspect differences among the samples, excluding
    non-contributing components can be helpful.

    Third, the exporter can ensure that transcripts belonging to the same ideal
    component are in a consistent orientation. When transcripts A and B belong
    to genes in the same ideal component, they can be in the same or reverse
    complement orientation. The exporter can reorient the transcript sequences
    so that all exported sequences from genes in the same ideal component are in
    the same orientation. In some cases, the exporter might not be able to
    reorient all sequences for transcripts belonging to genes in an ideal
    component so that they are all in the same orientation. This situation is
    caused by "inconsistencies" in the graph built to orient the
    transcripts. When an inconsistency is found, OrthologExporter raises an
    InconsistentGraphError unless the allow_inconsistent option is True. In that
    case, OrthologExporter will run a slower algorithm that attempts to best
    reorient the transcripts despite the inconsistency.

    Attributes:
        samples (list):                     List of samples in the analysis.
        ideal (list):                       List of ideal components.
        sample_gene_to_component (dict):    Mapping from (sample, gene) tuples
                                            to ideal component IDs.
        parse_transcript_id:                Function to parse FASTA sequence IDs
                                            into TranscriptID objects.
        ideal_ids (set):                    Ideal component IDs.
        strand_graph:                       Strand graph encoding relative
                                            orientations of transcripts.
        node_to_component_component (dict): Mapping from strand graph nodes to
                                            meta-strand graph components.
    """ 
    def __init__(
            self,
            sim: SampleSimilarity,
            parse_transcript_id: Callable[[str], TranscriptID],
            non_contributing: bool = True,
            consistent_strands: bool = True,
            allow_inconsistent: bool = True,
            jobs: int = 1,
            debug: bool = False,
    ):
        """Construct OrthologExporter for a SampleSimilarity with given options.

        This constructor accepts a few parameters that influence the behavior of
        the exporter.

        First, when non_contributing is False, ideal components in which all
        aligned regions are identical are not exported. These components are
        called "non-contributing components" because they introduce no
        differences.

        Second, when consistent_strands is True, the exporter will attempt to
        put all transcripts belonging to genes in the same ideal component in
        the same orientation.

        Putting all transcripts belonging to genes in a given ideal component
        can fail under certain circumstances. If the allow_inconsistent
        parameter is set to False, this constructor will fail in such a
        situtation, raising an InconsistentGraphError. If allow_inconsistent is
        instead set to True, the constructor will attempt to reorient the
        transcripts as well as it possibly can and continue despite the detected
        error(s).

        Parameters:
            sim:                       SampleSimilarity for the analysis.
            parse_transcript_id:       Function to parse transcript FASTA IDs.
            non_contributing (bool):   Admit non-contributing ideal components.
            consistent_strands (bool): Reorient transcripts per ideal component.
            allow_inconsistent (bool): Continue if strand graph is inconsistent.
            jobs (int):                Number of parallel jobs to use.
            debug (bool):              Enable debug behavior.
        """
        self.samples = sim.samples
        self.ideal = list(get_ideal_components(sim.graph, sim.sample_count))
        self.sample_gene_to_component = get_sample_gene_to_component(self.ideal)
        #print(self.sample_gene_to_component)
        self.parse_transcript_id = parse_transcript_id
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
                    self.parse_transcript_id,
                    jobs=jobs
                )
            dfs_assign_strands(self.strand_graph)
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
                msg = (
                    "Found {} strand mismatches in {} component groups.".format(
                        len(mismatches),
                        len(mismatch_component_components)
                    )
                )
                if not allow_inconsistent:
                    raise InconsistentGraphError(msg)
                # else:
                #     eprint("Attempting fix.")
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
                # eprint(
                #     "Reduced mismatches: {} -> {}".format(len(mismatches), mm2)
                # )
        #self.non_contributing = non_contributing
        #self.collapse = collapse
        #self.valid_tuples = set(map(tuple, sim.valid.itertuples(index=False)))
        #self.samples = sim.valid["sample"].drop_duplicates()

    def _name_ideal(self, sample: str, gene: int, isoform: int):
        """Return a string naming the ideal component a transcript belongs to.

        Although the isoform parameter is included to match the interface
        expected by other functions like concat_names, it is unused here.

        Parameters:
            sample (str):  Sample ID for the transcript
            gene (int):    Gene ID for the transcript
            isoform (int): (Unused) Isoform ID for the transcript

        Returns:
            A string representing the transcript's ideal component.
        """
        return "ideal_component_{}".format(
            self.sample_gene_to_component[(sample, gene)]
        )

    def _orient(
            self,
            t: tuple[str, int, int, Bio.SeqRecord]
    ) -> tuple[str, int, int, Bio.SeqRecord]:
        """Reorient the sequence according to the strand graph assignment.

        Parameters:
            t (tuple): Sample, gene, and isoform IDs and SeqRecord.

        Returns:
            Same sample, gene, and isoform IDs and reoriented SeqRecord.
        """
        try:
            if self.strand_graph.nodes[t[:-1]]["strand"] == -1:
                return t[:-1] + (
                    named_reverse_complement(
                        t[-1]
                    ),
                )                    
        except AttributeError:
            pass
        # except KeyError:
        #     from IPython import embed; embed()
        return t

    def by_sample(
            self,
            out_dir: Path,
            rename: Optional[Callable[[str, int, int], str]] = None,
            order: str = "after",
            make_all: bool = True
    ) -> dict[str, Path]:
        """Export orthologs, making one FASTA file per sample.

        Optionally, a function can be provided to rename the exported
        sequences. The function should accept the sample, gene, and isoform IDs
        of the original and return a string, the new FASTA ID/name.

        Parameters:
            out_dir:         Directory in which to create exported FASTA files.
            rename:          Optional function to rename sequences.
            order (str):     Put new name before or after the original.
            make_all (bool): Combine all files into all_ideal.fasta.

        Returns:
            A dictionary mapping sample IDs to paths to exported FASTA files.
        """
        if rename is None:
            rename = concat_names(self._name_ideal, order=order)
        sample_paths = {}
        for sample in self.samples:
            out_fn = out_dir / "{}_orthologs.fasta".format(
                path_to_sample(sample)
            )
            sample_paths[sample] = out_fn
            Bio.SeqIO.write(
                (
                    t[-1] for t in
                    renamed_seqs(
                        rename,
                        sorted(
                            (
                                self._orient(t)
                                for t in seq_tuples(
                                    sample,
                                    self.parse_transcript_id
                                )
                                if self.sample_gene_to_component.get(t[:-2])
                                in self.ideal_ids
                            ),
                            key=lambda x: self.sample_gene_to_component[x[:-2]]
                        ),
                    )
                ),
                out_fn,
                "fasta"
            )
            
        if make_all:
            self.make_all_ideal(sample_paths, out_dir)
        return sample_paths

    def by_component(
            self,
            out_dir: Path,
            rename: Optional[Callable[[str, int, int], str]] = None,
            order: str = "after",
            set_rlimit: bool = True,
            make_all: bool = True
    ) -> dict[int, Path]:
        """Export orthologs, making one FASTA file per ideal component.

        Optionally, a function can be provided to rename the exported
        sequences. The function should accept the sample, gene, and isoform IDs
        of the original and return a string, the new FASTA ID/name.

        Since there might be hundreds or even thousands of ideal components, it
        is easy for this function to exceed the soft open file resource limit
        imposed by the operating system. If the set_limit parameter is True,
        this method will try to automatically raise the open file limit to be
        above the expected number of open files. Since the new limit is only
        based on an estimate of the number of open file handles needed
        throughout this function, it is still possible for this method to exceed
        the rlimit. This is unlikely to happen for the
        rna_clique.export_orthologs command, but other programs calling this
        function should be careful to ensure the rlimit is sufficiently high.

        Parameters:
            out_dir:           Directory in which to save FASTA files.
            rename:            Optional function to rename sequences.
            order (str):       Put new name before or after the original.
            set_rlimit (bool): Try to increase rlimit.
            make_all (bool):   Combine all files into all_ideal.fasta.

        Returns:
            A dictionary mapping ideal component IDs exported FASTA file paths.
        """
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
                            t for t in seq_tuples(
                                sample,
                                self.parse_transcript_id
                            )
                            if t[:-2] in self.sample_gene_to_component
                        )
                ):
                    Bio.SeqIO.write(
                        self._orient((sample, gene, isoform, seq))[-1],
                        component_files[
                            self.sample_gene_to_component[
                                (
                                    sample,
                                    self.parse_transcript_id(seq.id).gene
                                )
                            ]
                        ],
                        "fasta"
                    )
        if make_all:
            self.make_all_ideal(component_paths, out_dir)
        return component_paths
        # if self.collapse:
        #     for f in component_paths.values():
        #         duplicates = defaultdict(list)
        #         for seq in Bio.SeqIO.parse(f, "r") as seq:
        #             duplicates[seq.seq].append(seq.description)

    def make_all_ideal(self, paths: Mapping[[Any], Path], export_out_dir: Path):
        """Create a file containing all exported transcripts.

        The created file is located at all_ideal.fasta under the provided
        export_out_dir directory and is a contains all sequences in the
        individual files located at the Paths given by the values of the paths
        argument.

        The combined file contains the same sequences as the individual files
        but changes their FASTA headers by appending the name of the file where
        each sequence was originally found.

        Parameters:
            paths (dict):   dict containing component FASTA files as values.
            export_out_dir: Directory in which to create combined file.

        """
        def seqs():            
            for path in paths.values():
                for seq in Bio.SeqIO.parse(path, "fasta"):
                    seq.description = ""
                    seq.title = ""
                    seq.name = ""
                    seq.id = "{}:{}".format(seq.id, path.stem)
                    #from IPython import embed; embed()
                    yield seq
        all_ideal_path = export_out_dir / "all_ideal.fasta"        
        Bio.SeqIO.write(seqs(), all_ideal_path, "fasta")
        #from IPython import embed; embed()        

def main():
    _, args, config = build_parser().get_arguments_and_config()
    config_module.RNACliqueConfigArgumentManager.make_output_dirs(config)
    args.export_output_dir.mkdir(exist_ok=True)
    sim = SampleSimilarity.from_filenames(
        config.graph,
        get_table_files(config.tables_dir),
        store_dfs=True
    )
    exporter = OrthologExporter(
        sim,
        TranscriptID.parser_from_re(config.transcript_id_regex),
        not args.remove_non_contributing,
        debug=args.debug,
        consistent_strands=not args.no_fix_strand,
        allow_inconsistent=args.allow_inconsistent,
        jobs=config.jobs,
    )
    getattr(exporter, "by_{}".format(args.by))(
        args.export_output_dir,
        order=args.concat_id_order,
        make_all=args.all
    )
    
if __name__ == "__main__":
    main()
