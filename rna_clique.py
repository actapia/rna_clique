import sys
import filtering_step
import multiprocessing
import filtered_distance
import itertools
from pathlib import Path

from transcripts import TranscriptID
from filtered_distance import SampleSimilarity
from similarity_computer import ComparisonSimilarityComputer
from transcripts import default_gene_re
from multiset_key_dict import MultisetKeyDict

out_dirs = filtering_step.out_dirs
out_files = filtering_step.out_files | {"output_matrix": "distance_matrix.h5"}
outputs = out_dirs | out_files

def build_parser():
    parser = filtering_step.build_parser()
    parser.add_argument(
        "--output-matrix",
        "-m",
        type=Path,
        help="Path to output distance matrix (overrides -O)."
    )
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

def handle_arguments():
    parser = build_parser()
    args = parser.parse_args()
    filtering_step.process_out_dir_args(parser, args, outputs)
    return args

def rna_clique(
        dirs,
        out_dir_1,
        out_dir_2,
        cache_dir,
        output_graph,
        output_matrix,
        top_genes,
        transcripts="transcripts.fasta",
        top_matches=1,
        id_parser=TranscriptID.parser_from_re(default_gene_re),
        evalue=1e-99,
        keep_all=True,
        store_dfs=False,
        jobs=multiprocessing.cpu_count() - 1,
):
    tables, table_paths, graph, num_tables = filtering_step.filtering_step(
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
    return sim


    
def main():
    args = handle_arguments()
    try:
        args.output_dir.mkdir(exist_ok=True)
    except AttributeError:
        pass
    for d in out_dirs:
        getattr(args, d).mkdir(exist_ok=True)
    id_parser = TranscriptID.parser_from_re(args.pattern)
    sim = rna_clique(
        args.dirs,
        args.out_dir_1,
        args.out_dir_2,
        args.cache_dir,
        args.output_graph,
        args.output_matrix,
        args.n,
        args.transcripts,
        args.N,
        id_parser,
        args.evalue,
        not args.no_keep_all,
        False,
        jobs=args.jobs
    )
    filtered_distance.writers[
        args.format
    ](
        sim.get_dissimilarity_df(),
        sys.stdout.buffer,
        header=args.header
    )
    
if __name__ == "__main__":
    main()
