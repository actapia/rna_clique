import shutil
import functools

import pysam
import Bio.Align
import networkx as nx

import config as config_module

from typing import Optional, Callable
from pathlib import Path
from collections import defaultdict, deque, namedtuple

from simple_blast import BlastDBCache, MultiformatBlastnSearch
from tqdm import tqdm

from gene_matches_tables import get_table_files
from filtered_distance import (
    SampleSimilarity,
    get_ideal_components,
)
from export_orthologs import build_strand_graph, get_sample_gene_to_component
from path_to_sample import path_to_sample
from transcripts import default_gene_re, TranscriptID

#default_gene_re = re.compile("^.*g([0-9]+)_i([0-9]+)")

def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager(
        description=(
            "BLAST search sequences of orthologous transcripts from ideal "
            "components."
        )
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
        help="Directory containing exported orthologs to search.",
    )
    arg_config.add_argument(
        "--all-ideal",
        "-a",
        type=Path,
        required=True,
        default={
            ("export_output_dir",): lambda export_output_dir:
            export_output_dir / "all_ideal.fasta"
        },
        help="FASTA file containing all sequences from ideal components.",
    )
    arg_config.add_argument(
        "--ortholog-db-cache",
        "-D",
        type=Path,
        default={
            ("export_output_dir",): lambda export_output_dir:
            export_output_dir / "db_cache"
        },
        help="Directory in which to store BLAST databases for orthologs.",
    )
    arg_config.add_argument(
        "--search-output-dir",
        "-S",
        type=Path,
        required=True,
        help="Output directory in which to store BLAST results.",
    )
    arg_config.add_argument(
        "--query",
        "-q",
        type=Path,
        required=True,
        help="FASTA file containing query sequences.",
    )
    arg_config.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug behavior.",
    )
    arg_config.add_argument(
        "--clean",
        action="store_true",
        help="Delete existing BLAST DB cache before beginning search.",
    )
    arg_config.add_argument(
        "--merge-sams",
        "-m",
        action="store_true",
        help="Merge extended search results into one file.",
    )
    arg_config.add_argument(
        "--extended-search",
        "-e",
        action="store_true",
        help="Search other isoforms of a gene that produces a hit.",
    )
    arg_config.add_argument(
        "--export-components",
        "-x",
        action="store_true",
        help="Save matching orientation graph components in extended search.",
    )
    return arg_config

SearchResult = namedtuple("SearchResult", ["hits", "seqs", "components"])

# TODO: Maybe move this to a class later to eliminate the need to construct a
# SampleSimilarity object.
def search(
        sim: SampleSimilarity,
        exported: Path,
        db_cache_loc: Path,
        out_dir: Path,
        query: Path,
        #gene_regex: re.Pattern = default_gene_re,
        parse_transcript_id: Callable[
            [str], TranscriptID
        ] = TranscriptID.parser_from_re(default_gene_re),
        path_to_sample: Optional[Callable[[Path], str]] = path_to_sample,
        # clean: bool = False,
        extended_search: bool = False,
        export_components: bool = True,
        merge_sams: bool = False,
        #strand_graph_out: tuple[nx.Graph, dict] = None
        strand_graph: nx.Graph = None,
        node_to_ccc: dict[tuple[str, int], nx.Graph] = None,
        jobs: int = 1,
        debug: bool = False,
) -> SearchResult:
    """Search for sequences within exported orthologs.

    This function requires a SampleSimilarity object from which it obtains the
    gene matches graph, gene matches tables, and sample count. The function also
    needs a Path to a file containing the concatenated ("all_ideal.fasta")
    exported orthologs, the output directory in which to store the search
    resutls, and a path to a FASTA file with the query sequences. The sequence
    IDs of the exported orthologs must be in a specific format---they should
    contain the original transcript FASTA ID followed by a colon (':') and then
    the sample ID. To use this format for an export (using OrthologExporter),
    use "after" for the order parameter to the by_components method.

    By default, the function searches the exported orthologs for the given query
    sequences, producing a SAM alignment named "queries.sam" in the provided
    output directory. The function also produces a "subjects.fasta" file
    containing the sequences of the transcripts for which BLAST obtained
    alignments with the query sequences.

    Optionally, an "extended search" can be performed by passing True for the
    extended_search parameter. In the extended search, when a transcript matches
    one of the query sequences in the initial search, additional searches are
    performed to compare the query sequences with other isoforms of the same
    gene. These additional searches are performed with a lower e-value threshold
    to possibly find hits missing in the original search. For each transcript
    appearing in the initial search, the additional isoforms searched are
    exactly those isoforms of the same gene that appear in the same strand graph
    connected component are search. Note that since two isoforms of the same
    gene might not appear in the same strand graph connected component, the
    presence of one isoform in the initial search does not guarantee that any
    other isoform will be searched in the extended search.

    Additionally, during the extended search, this function can export the
    components of the strand graph corresponding to the ideal components to
    which transcripts appearing in the BLAST results belong. This behavior can
    be enabled by passing a value of True for the export_component parameter.

    When this function searches additional isoforms in the extended search,
    every isoform's search results are written to separate SAM file with the
    pattern "{sample_id}_g{gene_id}_i{isoform_id}.sam". Sometimes it is more
    convenient to have all of these results in a single SAM file. To enable
    creation of such a "graph.sam" file, pass True for the merge_sams parameter.

    One goal of the extended search is to help the user determine if a sequence
    found among the exported orthologs is really present or if it might be an
    assembly artifact. In the latter case, one might see that matches occur only
    in one isoform, and such an isoform might be disconnected from the others in
    the strand graph if it is of sufficiently poor quality.

    In addition to the options above related to the extended search, this
    function optionally accepts a function for parsing transcript IDs.

    By default, this function will construct a strand graph and a dict mapping
    ideal component nodes to their corresponding components in the meta-strand
    graph. Since the OrthologExporter computes these values, it can be more
    efficient to reuse them when both exporting and searching are being
    performed as part of the same code. This function accepts the strand graph
    and dict via the strand_graph and node_to_ccc parameters, respectively.

    This function returns a SearchResult object, which is a namedtuple with
    three attributes: hits, seqs, and components. The first, hits, is the number
    of hits found in the initial search. The second is the number of matching
    sequences in the initial search. The third and last is the number of ideal
    components (or, equivalently, meta-strand graph components) in which there
    were matching sequences.

    Parameters:
        sim:                      SampleSimilarity for the RNA-clique analysis.
        exported:                 Path to concatenated exported orthologs.
        db_cache_loc:             Intermediate directory for BLAST DBs.
        out_dir:                  Output directory for search results.
        query:                    Path to sequences to search for.
        parse_transcript_id:      FASTA ID to TranscriptID parsing function.
        path_to_sample:           Top genes path to sample name function.
        extended_search (bool):   Search other isotigs of matched genes.
        export_components (bool): Save matching strand graph components in
                                  extended search.
        merge_sams (bool):        Merge extended search results into one SAM.
        strand_graph:             Strand graph for the analysis.
        node_to_ccc (dict):       dict mapping ideal component nodes to
                                  meta-strand graph connected components.
        jobs (int):               Number of parallel jobs to use.
        debug (bool):             Enable debug behavior.

    Returns:
        The SearchResult object representing results of searching the exports.
    """
    db_cache_loc.mkdir(exist_ok=True)
    out_dir.mkdir(exist_ok=True)
    # TODO: Why use a DB cache here?
    cache = BlastDBCache(db_cache_loc)
    exports = [exported]
    cache.makedb(exports)
    # Start by just searching for the query sequences in all of the exported
    # transcripts.
    search = MultiformatBlastnSearch(query, exports, db_cache=cache)
    tab_search = search.to_search(6)
    if not tab_search.hits.empty:
        # print("Found {} hits.".format(tab_search.hits.shape[0]))
        sample_to_path = {path_to_sample(v): v for v in sim.samples}
        # Save the initial search alignments in SAM format.
        Bio.Align.write(
            search.to_sam(subject_as_reference=True).hits,
            out_dir / "queries.sam",
            "sam"
        )
        subjects = set()
        export_index = Bio.SeqIO.index(exported, "fasta")
        ideal = list(get_ideal_components(sim.graph, sim.sample_count))
        sample_gene_to_component = get_sample_gene_to_component(ideal)
        # TODO: See if we can avoid rebuilding node_to_ccc when only
        # strand_graph is provided.
        if strand_graph is None or node_to_ccc is None:
            strand_graph, node_to_ccc = build_strand_graph(
                sim,
                sample_gene_to_component,
                parse_transcript_id,
                jobs=jobs
            )
        # Map meta-strand graph components to lists of nodes they contain.
        # cccs ONLY contains the components for which the corresponding
        # nodes can be found in the BLAST results.
        cccs = defaultdict(list)
        for full_seq_id in tab_search.hits["sseqid"].drop_duplicates():
            seq_id, sample = full_seq_id.split(":")
            node = (sample_to_path[sample],) + \
                parse_transcript_id(seq_id)[1:]
            cccs[node_to_ccc[node]].append(node)
            subjects.add(full_seq_id)
        sam_paths = []
        if extended_search:
            # Get a mapping from gene matches graph nodes to sequence IDs.
            node_to_seq_id = {}
            for full_seq_id in export_index:
                seq_id, sample = full_seq_id.split(":")
                node = (sample_to_path[sample],) + \
                    parse_transcript_id(full_seq_id)[1:]
                node_to_seq_id[node] = full_seq_id
            # Map sample-gene pairs to ideal component indices.
            # print("Going over component connected components.")
            for ccc, nodes in tqdm(cccs.items()):
                if export_components:
                    # Get the corresponding ideal component index for the
                    # nodes. All of them are guarnteed to be in the same ideal
                    # component because they all come from the same meta-strand
                    # graph component.
                    ideal_index = sample_gene_to_component[
                        next(
                            iter(nodes)
                        )[:-1]
                    ]
                    # Write the strand graph subgraph corresponding to the ideal
                    # component.
                    nx.write_graphml(
                        functools.reduce(nx.union, ccc.nodes),
                        out_dir / f"ideal_component_{ideal_index}.graphml"
                    )
                for node in nodes:
                    # Get the strand graph connected component corresponding to
                    # the node. A more efficient way of doing this might be
                    # needed eventually.
                    cc = next(x for x in ccc.nodes if node in x.nodes)
                    # Perform a search in the strand graph for nodes from the
                    # same gene and perform additional searches with lower
                    # e-value thresholds for those nodes.
                    seen = {node}
                    to_search = deque()
                    to_search.append((None, node))
                    while to_search:
                        prev, n = to_search.pop()
                        subjects.add(node_to_seq_id[n])
                        same_neighbors = {
                            m for m in cc.neighbors(n) if m[0] == node[0]
                        }
                        to_search.extend(
                            (n, m) for m in same_neighbors
                            if m not in seen
                        )
                        seen |= same_neighbors
                        queries = [
                            export_index[node_to_seq_id[m]]
                            for m in cc.neighbors(n)
                            if m != prev
                        ]
                        if queries:
                            with MultiformatBlastnSearch.from_sequences(
                                    [export_index[node_to_seq_id[n]]],
                                    queries,
                                    evalue=1e-5,
                            ) as new_search:
                                out_path = out_dir / "{}_g{}_i{}.sam".format(
                                    path_to_sample(n[0]),
                                    *n[1:]
                                )
                                Bio.Align.write(
                                    new_search.to_sam(
                                        subject_as_reference=False
                                    ).hits,
                                    out_path,
                                    "sam"
                                )
                                sam_paths.append(out_path)
        if sam_paths and merge_sams:
            pysam.samtools.merge(
                "-o",
                str(out_dir / "graph.sam"),
                *map(str, sam_paths)
            )
        Bio.SeqIO.write(
            (
                export_index[x]
                for x in subjects
            ),
            out_dir / "subjects.fasta",
            "fasta"
        )
        return_result = SearchResult(
            tab_search.hits.shape[0],
            len(tab_search.hits["sseqid"].drop_duplicates()),
            len(cccs)
        )
        # print("Returning", return_result)
        return return_result
    else:
        return SearchResult(0, 0, 0)

def main():
    _, args, config = build_parser().get_arguments_and_config()
    args.ortholog_db_cache.mkdir(exist_ok=True)
    args.search_output_dir.mkdir(exist_ok=True)
    if args.ortholog_db_cache.exists() and args.clean:
        shutil.rmtree(args.ortholog_db_cache)
    sim = SampleSimilarity.from_filenames(
        config.graph,
        get_table_files(config.tables_dir),
        store_dfs=True,
    )
    search(
        sim,
        exported=args.all_ideal,
        db_cache_loc=args.ortholog_db_cache,
        out_dir=args.search_output_dir,
        query=args.query,
        parse_transcript_id=TranscriptID.parser_from_re(
            config.transcript_id_regex
        ),
        extended_search=args.extended_search,
        export_components=args.export_components,
        merge_sams=args.merge_sams,
        jobs=config.jobs,
        debug=args.debug,
    )

if __name__ == "__main__":
    main()
