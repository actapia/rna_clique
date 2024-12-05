import argparse
import re
import pickle
import pandas as pd
import itertools

from pathlib import Path

from tqdm import tqdm

from subset_comparisons import (
    handle_filters,
    matcher,
    make_subset_comparisons
)
from path_to_sample import sample_re
from build_graph import build_graph
from gene_matches_tables import read_table

from collections.abc import Iterable
from typing import Iterator



def handle_arguments():
    parser = argparse.ArgumentParser(
        description="subset a previous analysis, creating a new graph"
    )
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
        "--exclude",
        "-x",
        nargs="+",
        default=[],
        help="samples to exclude (default is none)"
    )
    parser.add_argument(
        "--include",
        "-y",
        nargs="+",
        default=[],
        help="samples to include (default is all)"
    )
    parser.add_argument(
        "--include-regex",
        "-Y",
        type=re.compile,
        help="regular expression specifying which sample names to include"
    )
    parser.add_argument(
        "--sample-name-regex",
        "-r",
        help="regular expression to apply to sample names",
        type=re.compile,
        default=sample_re
    )
    parser.add_argument(
        "--include-file",
        type=Path,
        help="file containing samples to include"
    )
    parser.add_argument(
        "--exclude-file",
        type=Path,
        help="file containing samples to exclude"
    )
    parser.add_argument(
        "--show-included",
        action="store_true",
        help="show which samples would be included and exit"
    )
    parser.add_argument(
        "--show-parsed-paths",
        action="store_true",
        help="show parsed paths"
    )
    return parser.parse_args()

def multi_glob(path: Path, globs: Iterable[str]) -> Iterator[Path]:
    return itertools.chain(*map(path.glob, globs))

def main():
    args = handle_arguments()
    include = handle_filters(args.include, args.include_file)
    if not include:
        include = None
    exclude = handle_filters(args.exclude, args.exclude_file)
    if not exclude:
        exclude = None
    matches = matcher(
        include,
        exclude,
        args.include_regex
    )
    inputs = list(multi_glob(args.input_dir / "od2", ["*.pkl", "*.h5"]))    
    if args.show_included or args.show_parsed_paths:
        for df_path in inputs:
            df = read_table(df_path, head=1)
            if args.show_parsed_paths:
                print(
                    *(
                        args.sample_name_regex.search(
                            Path(df[x + "sample"].iloc[0]).name
                        ).group(1)
                        for x in ["q", "s"]                
                    )
                )
            if args.show_included and all(
                matches(
                    args.sample_name_regex.search(
                        Path(df[x + "sample"].iloc[0]).name
                    ).group(1)
                )
                for x in ["q", "s"]
            ):
                print(df_path)

    else:
        args.output_dir.mkdir(exist_ok=True)
        od2 = args.output_dir / "od2"
        od2.mkdir(exist_ok=True)

        graph = build_graph(
            make_subset_comparisons(
                tqdm(inputs),
                od2,
                matches,
                args.sample_name_regex
            )
        )
        with open(args.output_dir / "graph.pkl", "wb") as f:
            pickle.dump(graph, f, pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    main()


    
