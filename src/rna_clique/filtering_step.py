import multiprocessing
import itertools
import pickle

from typing import Iterable

import networkx as nx
import pandas as pd

from pathlib import Path
from typing import Callable

from joblib import Parallel, delayed
from tqdm import tqdm

from . import app
from . import config as config_module
from .transcripts import default_gene_re, TranscriptID, TranscriptIDParseError
from .select_top_genes_all import select_top_and_save
from .find_all_pairs import find_all_pairs
from .build_graph import build_graph
from .similarity_computer import ComparisonSimilarityComputer
from .app import set_except_hook, validate_input_dirs

def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager(
        description="Select top genes and get gene matches tables and graph."
    )
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


# TODO: Link to a better explanation of the parameters rather than trying to
# rehash the whole algorithm here.
def filtering_step(
        dirs: Iterable[Path],
        out_dir_1: Path,
        out_dir_2: Path,
        cache_dir: Path,
        output_graph: Path,
        top_genes: int,
        transcripts: str = "transcripts.fasta",
        top_matches: int = 1,
        id_parser: Callable[
            [str],
            TranscriptID
        ]=TranscriptID.parser_from_re(default_gene_re),
        evalue: float = 1e-99,
        keep_all: bool = True,
        jobs: int = multiprocessing.cpu_count() - 1,
) -> tuple[Iterable[pd.DataFrame], Iterable[Path], nx.Graph]:
    """Perform the filtering step (phase 1) of RNA-clique.

    This function performs the full filtering step of RNA-clique, which has also
    been referred to as "phase 1" in some RNA-clique documentation. Phase 1
    consists of the following steps:

        1. Selection of top n genes by k-mer coverage for every sample.
        2. BLASTn searching for top n genes of every ordered sample pair.
        3. Processing of the BLASTn searches to get gene matches tables.
        4. Building the gene matches graph.

    The third step involves finding for each unordered pair of samples the best
    matches in both directions for every gene belonging to one of the samples. A
    pair of genes, one from each sample, makes it into the "gene matches table"
    for the pair of samples only when that pair of genes appears among the top N
    matches for the first sample in both directions. This parameter is usually
    set to 1, so a match is only counted if it is (one of) the best in both
    directions. In the rare case of ties, you can keep all mathces using the
    keep_all parameter.

    For a more thorough explanation, please refer to the RNA-clique paper.

    This function mainly performs I/O, but it also returns three objects that
    are convenient for downstream processing. First, the function returns an
    iterable of the gene matches tables. Second, the function returns an
    iterable of paths to the gene matches tables. Third, the function returns
    the gene matches graph.

    To reduce memory requirements, the function avoids loading all gene matches
    tables into memory at once. To this end, the returned iterable over gene
    matches tables does not iterate over tables stored in memory. Instead,
    tables are loaded from disk as the iteratable is iterated.

    Parameters:
        dirs:              Input directory containing transcriptomes.
        out_dir_1:         Output directory for storing top genes by coverage.
        out_dir_2:         Output directory for storing gene matches tables.
        cache_dir:         Intermediate directory storing BLAST DB caches.
        output_graph:      Path to output gene matches graph pickle.
        top_genes (int):   Number of top genes to select.
        transcripts (str): Name of transcript FASTA files within input dirs.
        top_matches (int): Threshold for counting matches between directions.
        id_parser:         Function to parse TranscriptIDs from FASTA IDs.
        evalue (float):    BLAST search e-value threshold.
        keep_all (bool):   Whether to keep all matches in case of a tie.
        jobs (int):        Number of parallel jobs to use.

    Returns:
        Two iterables with gene matches tables and paths, gene matches graph.
    """
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
    with set_except_hook():
        _, args, config = build_parser().get_arguments_and_config()
    with set_except_hook(config.verbose):
        config_module.RNACliqueConfigArgumentManager.make_output_dirs(config)
        id_parser = TranscriptID.parser_from_re(config.transcript_id_regex)
        validate_input_dirs(config)
        try:
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
        except TranscriptIDParseError:
            app.print_transcript_id_parse_error_message(
                config.transcript_id_regex
            )
            raise        
        config.path_to_sample = pts
        config.mark_finish()
        if args.output_config:
            config.yaml_save(args.output_config)

if __name__ == "__main__":
    main()
