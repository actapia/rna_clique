import argparse
import re

from pathlib import Path

import pandas as pd

from tqdm import tqdm

from build_graph import build_graph

default_filter_regex = re.compile("(.*)")

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_arguments(
        "-O",
        "--output-dir",
        type=Path,
        required=True,
        help="path to the directory in which to store filtered comparisons"
    )
    parser.add_argument(
        "-i",
        "--inputs",
        nargs="+",
        type=Path,
        required=True,
        help="pickles containing the gene matches tables"
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

def matcher(included, filter_regex):
    def inner(x):
        return x in included or (filter_regex and filter_regex.search(x))
    return inner

def main():
    args = handle_arguments()
    args.out_dir.mkdir(exist_ok=True)
    include = set(args.filter)
    try:
        with open(args.filter_file, "r") as filter_file:
            include |= set(l.rstrip() for l in filter_file)
    except TypeError:
        pass
    matches = matcher(include, args.filter_regex)
    for df_path in tqdm(args.inputs):
        df = pd.read_pickle(df_path)
        if all(
                matches(
                    args.file_regex.search(
                        df[x + "sample"].iloc[0]
                    ).group(1)
                    for x in ["q", "s"]
                )
        ):
            dest = args.out_dir / df_path.name
            dest.symlink_to(df_path.relative_to(dest))
            
if __name__ == "__main__":
    main()
