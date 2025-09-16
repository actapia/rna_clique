import multiprocessing
import sys

from pathlib import Path
from typing import Callable, Iterable

from multiset_key_dict import MultisetKeyDict

from .transcripts import TranscriptID
from .filtered_distance import SampleSimilarity
from .similarity_computer import ComparisonSimilarityComputer
from .transcripts import default_gene_re
from . import filtered_distance, filtering_step
from . import config as config_module

def build_parser():
    parser = filtering_step.build_parser()
    parser.expose_fields_with_default_aliases("matrix")
    parser.add_argument(
        "-f",
        "--format",
        choices=filtered_distance.writers,
        default="table",
        help="Format for writing distance matrix to stdout."
    )
    parser.add_argument(
        "--header",
        action="store_true",
        help="Include header in distance matrix written to stdout."
    )
    return parser

def rna_clique(
        dirs: Iterable[Path],
        out_dir_1: Path,
        out_dir_2: Path,
        cache_dir: Path,
        output_graph: Path,
        output_matrix: Path,
        top_genes: int,
        transcripts: str = "transcripts.fasta",
        top_matches: int = 1,
        id_parser: Callable[
            [str],
            TranscriptID
        ] = TranscriptID.parser_from_re(default_gene_re),
        evalue: float = 1e-99,
        keep_all: bool = True,
        store_dfs: bool = False,
        jobs: int = multiprocessing.cpu_count() - 1,
) -> tuple[SampleSimilarity, dict[Path, str]]:
    tables, table_paths, graph, num_tables, pts = filtering_step.filtering_step(
        dirs,
        out_dir_1,
        out_dir_2,
        cache_dir,
        output_graph,
        top_genes,
        transcripts,
        top_matches,
        id_parser,
        evalue,
        keep_all,
        jobs,
    )
    tables = ComparisonSimilarityComputer.mapping_from_dfs(tables)
    if store_dfs:
        tables = MultisetKeyDict(tables)
    sim = SampleSimilarity(
        graph,
        tables,
    )
    mat = sim.get_dissimilarity_df()
    mat.to_hdf(output_matrix, "matrix")
    return sim, pts
    
def main():
    _, args, config = build_parser().get_arguments_and_config()
    config_module.RNACliqueConfigArgumentManager.make_output_dirs(config)
    id_parser = TranscriptID.parser_from_re(config.transcript_id_regex)
    sim, pts = rna_clique(
        config.input_dirs,
        config.top_genes_dir,
        config.tables_dir,
        config.cache_dir,
        config.graph,
        config.matrix,
        config.top_genes,
        config.transcripts_name,
        config.top_matches,
        id_parser,
        config.evalue,
        config.keep_all,
        False,
        jobs=config.jobs
    )
    filtered_distance.writers[
        args.format
    ](
        sim.get_dissimilarity_df(),
        sys.stdout.buffer,
        header=args.header
    )
    config.path_to_sample = pts
    config.mark_finish()
    config.yaml_save(args.output_config)

    
if __name__ == "__main__":
    main()
