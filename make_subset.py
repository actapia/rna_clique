import argparse
import re
import pickle

from pathlib import Path

from tqdm import tqdm

from subset_comparisons import (
    default_filter_regex,
    handle_filters,
    matcher,
    make_subset_comparisons
)
from build_graph import build_graph

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-O",
        "--output-dir",
        type=Path,
        required=True,
        help="directory in which to store subset comparisons and graph"
    )
    parser.add_argument(
        "-I",
        "--input-dir",
        type=Path,
        required=True,
        help="directory for full set of data"
    )
    parser.add_argument(
        "--filter",
        "-f",
        nargs="+",
        default=[],
        help="samples to include (if not provided, the problem includes all)"
    )
    parser.add_argument(
        "--filter-regex",
        "-R",
        type=re.compile,
        help="regular expression specifying which sample names to include"
    )
    parser.add_argument(
        "--sample-name-regex",
        "-r",
        help="regular expression to apply to sample names",
        type=re.compile,
        default=default_filter_regex
    )
    parser.add_argument(
        "--filter-file",
        "-F",
        type=Path,
        help="file containing sample to include"
    )
    return parser.parse_args()

def main():
    args = handle_arguments()
    args.output_dir.mkdir(exist_ok=True)
    od2 = args.output_dir / "od2"
    od2.mkdir(exist_ok=True)
    matches = matcher(
        handle_filters(args.filter, args.filter_file),
        args.filter_regex
    )
    graph = build_graph(
        make_subset_comparisons(
            tqdm(list((args.input_dir / "od2").glob("*pkl"))),
            od2,
            matches,
            args.sample_name_regex
        )
    )
    with open(args.output_dir / "graph.pkl", "wb") as f:
        pickle.dump(graph, f, pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    main()


    
