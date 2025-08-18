import argparse
import multiprocessing
import re
import itertools
import pickle

from typing import Iterable

import networkx as nx
import pandas as pd
import Bio.SeqIO

from pathlib import Path
from joblib import Parallel, delayed
from tqdm import tqdm

from transcripts import default_gene_re, TranscriptID
from select_top_genes import TopGeneSelector
from find_all_pairs import find_all_pairs
from build_graph import build_graph
from similarity_computer import ComparisonSimilarityComputer

out_dirs = {
    "out_dir_1": "od1",
    "out_dir_2": "od2",
    "cache_dir": "db_cache",
}
out_files = {
    "output_graph": "graph.pkl"
}
outputs = out_dirs | out_files

def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--top-genes",
        "-n",
        type=int,
        help="Top n genes to select by k-mer coverage.",
        dest="n",
        required=True
    )
    parser.add_argument(
        "--top-matches",
        "-N",
        type=int,
        help="Count a match if it is within top N in both directions.",
        dest="N",
        default=1
    )
    parser.add_argument(
        "--transcripts",
        "-t",
        help="Name of transcripts files.",
        default="transcripts.fasta"
    )
    parser.add_argument(
        "--out-dir-1",
        "-O1",
        type=Path,
        help="Intermediate out directory containing top genes (overrides -O)."
    )
    parser.add_argument(
        "--out-dir-2",
        "-O2",
        type=Path,
        help="Intermediate out directory containing matches (overrides -O)."
    )
    parser.add_argument(
        "--pattern",
        "-p",
        type=re.compile,
        help="Regular expression for parsing transcript FASTA headers.",
        default=default_gene_re,
    )
    parser.add_argument(
        "--evalue",
        "-e",
        type=float,
        help="Cutoff evalue to use in BLAST searches.",
        default=1e-99
    )
    parser.add_argument(
        "--no-keep-all",
        action="store_true",
        help="Do not keep all matches in case of a tie."
    )
    parser.add_argument(
        "--output-graph",
        "-g",
        type=Path,
        help="Path to output graph (overrides -O)."
    )
    parser.add_argument(
        "--jobs",
        "-j",
        type=int,
        help="Number of parallel jobs to use.",
        default=multiprocessing.cpu_count() - 1
    )
    parser.add_argument(
        "--cache-dir",
        "-C",
        type=Path,
        help="Directory in which to store BLAST DBs (overrides -O)."
    )
    parser.add_argument(
        "--output-dir",
        "-O",
        type=Path,
        help="Output directory."
    )
    parser.add_argument(
        "dirs",
        type=Path,
        nargs="+"
    )
    return parser

def handle_arguments():
    parser = build_parser()
    args = parser.parse_args()
    process_out_dir_args(parser, args, outputs)
    return args

def process_out_dir_args(parser, args, outs):
    for arg, f in outs.items():
        if getattr(args, arg, None) is None:
            setattr(args, arg, args.output_dir / f)
    try:
        missing = next(x for x in outs if getattr(args, x) is None)
        parser.error(
            "Must provide --output-dir or --{}".format(missing.replace("_","-"))
        )
    except StopIteration:
        pass


def select_top_and_save(out_dir, transcripts, x: Path, *args):
    out = out_dir / (x.stem + "_top.fasta")
    Bio.SeqIO.write(
        TopGeneSelector.from_path(
            x / transcripts,
            *args
        ).get_top_gene_seqs(),
        out,
        "fasta"
    )
    return (out, x.stem)

def filtering_step(
        dirs,
        out_dir_1,
        out_dir_2,
        cache_dir,
        output_graph,
        top_genes,
        transcripts="transcripts.fasta",
        top_matches=1,
        id_parser=TranscriptID.parser_from_re(default_gene_re),
        evalue=1e-99,
        keep_all=True,
        jobs=multiprocessing.cpu_count() - 1,
) -> tuple[Iterable[pd.DataFrame], Iterable[Path], nx.Graph]:
    path_to_sample = dict(
        Parallel(n_jobs=jobs)(
            map(
                delayed(
                    lambda x: select_top_and_save(
                        out_dir_1,
                        transcripts,
                        x,
                        top_genes,
                        id_parser
                    )
                ),
                dirs
            )
        )
    )
    tables, table_paths, num_tables = find_all_pairs(
        path_to_sample,
        out_dir_2,
        cache_dir,
        path_to_sample.__getitem__,
        hf_args=[
            id_parser,
            top_matches,
            evalue,
            keep_all
        ],
        jobs=jobs,
    )
    graph = build_graph(tqdm(tables, total=num_tables))
    with open(output_graph, "wb") as f:
        pickle.dump(graph, f, pickle.HIGHEST_PROTOCOL)
    table_paths1, table_paths2 = itertools.tee(table_paths)
    return map(
        ComparisonSimilarityComputer._read_table,
        table_paths1
    ), table_paths2, graph, num_tables
    

def main():
    args = handle_arguments()
    try:
        args.output_dir.mkdir(exist_ok=True)
    except AttributeError:
        pass
    for d in out_dirs:
        getattr(args, d).mkdir(exist_ok=True)
    id_parser = TranscriptID.parser_from_re(args.pattern)
    filtering_step(
        args.dirs,
        args.out_dir_1,
        args.out_dir_2,
        args.cache_dir,
        args.output_graph,
        args.n,
        args.transcripts,        
        args.N,
        id_parser,
        args.evalue,
        not args.no_keep_all,
        args.jobs
    )


if __name__ == "__main__":
    main()
