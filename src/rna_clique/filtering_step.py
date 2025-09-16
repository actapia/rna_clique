import multiprocessing
import itertools
import pickle

import networkx as nx
import pandas as pd

from pathlib import Path
from typing import Iterable

from joblib import Parallel, delayed
from tqdm import tqdm

from .transcripts import default_gene_re, TranscriptID
from .select_top_genes_all import select_top_and_save
from .find_all_pairs import find_all_pairs
from .build_graph import build_graph
from .similarity_computer import ComparisonSimilarityComputer
from . import config as config_module

def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager()
    arg_config.expose_fields_with_default_aliases(
        "top_genes",
        "top_matches",
        "transcripts_name",
        "top_genes_dir",
        "tables_dir",
        "transcript_id_regex",
        "evalue",
        "jobs",
        "cache_dir",
        "graph",
        required=True
    )
    arg_config.expose_fields_with_default_aliases("output_dir", "title")
    arg_config.add_argument(
        "--no-keep-all",
        dest="keep_all",
        action="store_false",
        help="Do not keep all matches in case of a tie."
    )
    arg_config.add_output_config_argument()
    arg_config.expose_config_field(
        "input_dirs",
        nargs="*",
        #const=None,
        positional=True,
    )
    return arg_config


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
    ), table_paths2, graph, num_tables, path_to_sample
    

def main():
    _, args, config = build_parser().get_arguments_and_config()
    config_module.RNACliqueConfigArgumentManager.make_output_dirs(config)
    id_parser = TranscriptID.parser_from_re(config.transcript_id_regex)
    pts = filtering_step(
        config.input_dirs,
        config.top_genes_dir,
        config.tables_dir,
        config.cache_dir,
        config.graph,
        config.top_genes,
        config.transcripts_name,
        config.top_matches,
        id_parser,
        config.evalue,
        config.keep_all,
        config.jobs
    )[-1]
    config.path_to_sample = pts
    config.mark_finish()
    config.yaml_save(args.output_config)

if __name__ == "__main__":
    main()
