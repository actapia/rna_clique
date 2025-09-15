import argparse
import re
import pysam
import shutil
import functools
import Bio.Align
import networkx as nx
import config as config_module
from collections.abc import Iterable
from typing import Optional, Callable
from pathlib import Path
from collections import defaultdict, deque, namedtuple
from gene_matches_tables import get_table_files

from simple_blast import BlastDBCache, MultiformatBlastnSearch

from filtered_distance import (
    SampleSimilarity,
    get_ideal_components,
)
from export_orthologs import build_strand_graph, get_sample_gene_to_component
from path_to_sample import path_to_sample
from transcripts import default_gene_re, TranscriptID

from tqdm import tqdm

#default_gene_re = re.compile("^.*g([0-9]+)_i([0-9]+)")

def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager()
    arg_config.expose_fields_with_default_aliases(
        "graph",
        "tables_dir",
        "jobs",
        "transcript_id_regex",
        required=True
    )
    arg_config.expose_config_field(
        "output_dir",
        aliases=["--analysis-root", "--rna-clique-output-dir", "-A"]        
    )
    arg_config.add_argument(
        "--export-output-dir",
        "-X",
        type=Path,
    )
    arg_config.add_argument(
        "--all-ideal",
        "-a",
        type=Path,
        required=True,
        default={
            ("export_output_dir",): lambda export_output_dir:
            export_output_dir / "all_ideal.fasta"
        }
    )

    arg_config.add_argument(
        "--ortholog-db-cache",        
        "-D",
        type=Path,
        default={
            ("export_output_dir",): lambda export_output_dir:
            export_output_dir / "db_cache"
        }
    )
    arg_config.add_argument(
        "--search-output-dir",
        "-S",
        type=Path,
        required=True
    )
    arg_config.add_argument(
        "--query",
        "-q",
        type=Path,
        required=True
    )
    arg_config.add_argument(
        "--debug",
        action="store_true"
    )
    arg_config.add_argument(
        "--clean",
        action="store_true"
    )
    arg_config.add_argument(
        "--merge-sams",
        "-m",
        action="store_true"
    )
    return arg_config

SearchResult = namedtuple("SearchResult", ["hits", "seqs", "components"])

def search(
        graph_loc: Path,
        comparisons_loc: Iterable[Path],
        exported: Path,
        db_cache_loc: Path,
        out_dir: Path,
        query: Path,
        debug: bool = False,
        sample_count: Optional[int] = None,
        #gene_regex: re.Pattern = default_gene_re,
        parse_transcript_id: Callable[
            [str], TranscriptID
        ] = TranscriptID.parser_from_re(default_gene_re),
        clean: bool = False,
        merge_sams: bool = False,
        jobs: int = 1,
        sim: SampleSimilarity = None,
        strand_graph_out: tuple[nx.Graph, dict] = None
):
    if db_cache_loc.exists() and clean:
        shutil.rmtree(db_cache_loc)
    db_cache_loc.mkdir(exist_ok=True)
    out_dir.mkdir(exist_ok=True)
    cache = BlastDBCache(db_cache_loc)
    exports = [exported]
    cache.makedb(exports)
    search = MultiformatBlastnSearch(query, exports, db_cache=cache)
    tab_search = search.to_search(6)
    if not tab_search.hits.empty:
        print("Found {} hits.".format(tab_search.hits.shape[0]))
        if sim is None:
            sim = SampleSimilarity.from_filenames(
                graph_loc,
                tqdm(comparisons_loc),
                sample_count=sample_count,
                store_dfs=True
            )
        sample_to_path = {path_to_sample(v): v for v in sim.samples}
        Bio.Align.write(
            search.to_sam(subject_as_reference=True).hits,
            out_dir / "queries.sam",
            "sam"
        )
        export_index = Bio.SeqIO.index(exported, "fasta")
        node_to_seq_id = {}
        for full_seq_id in export_index:
            seq_id, sample = full_seq_id.split(":")
            node = (sample_to_path[sample],) + \
                parse_transcript_id(full_seq_id)[1:]
            node_to_seq_id[node] = full_seq_id
        ideal = list(get_ideal_components(sim.graph, sim.sample_count))
        sample_gene_to_component = get_sample_gene_to_component(ideal)
        if strand_graph_out is None:
            strand_graph, node_to_ccc = build_strand_graph(
                sim,
                sample_gene_to_component,
                parse_transcript_id,
                jobs=jobs
            )
        else:
            strand_graph, node_to_ccc = strand_graph_out
        cccs = defaultdict(list)
        for seq_id in tab_search.hits["sseqid"].drop_duplicates():
            seq_id, sample = seq_id.split(":")
            node = (sample_to_path[sample],) + \
                parse_transcript_id(seq_id)[1:]
            cccs[node_to_ccc[node]].append(node)
        print("Going over component connected components.")
        subjects = set()
        sam_paths = []
        for ccc, nodes in tqdm(cccs.items()):
            ideal_index = sample_gene_to_component[next(iter(nodes))[:-1]]
            nx.write_graphml(
                functools.reduce(nx.union, ccc.nodes),
                out_dir / "ideal_component_{}.graphml".format(ideal_index)
            )
            for node in nodes:
                cc = next(x for x in ccc.nodes if node in x.nodes)
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
        if merge_sams:
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
        print("Returning", return_result)
        return return_result
    else:
        return SearchResult(0, 0, 0)

def main():
    _, args, config = build_parser().get_arguments_and_config()
    args.ortholog_db_cache.mkdir(exist_ok=True)
    args.search_output_dir.mkdir(exist_ok=True)
    search(
        graph_loc=config.graph,
        comparisons_loc=get_table_files(config.tables_dir),
        exported=args.all_ideal,
        db_cache_loc=args.ortholog_db_cache,
        out_dir=args.search_output_dir,
        query=args.query,
        debug=args.debug,
        #gene_regex=TranscriptID.parser_from_re(config.transcript_id_regex),
        clean=args.clean,
        merge_sams=args.merge_sams,
        jobs=config.jobs
    )
        
if __name__ == "__main__":
    main()
