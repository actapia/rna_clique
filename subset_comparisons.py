import argparse
import os
import re

from pathlib import Path

import pandas as pd

from tqdm import tqdm

from build_graph import build_graph

from more_itertools import consume

default_filter_regex = re.compile("(.*)")

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
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

def relative_to(p1: Path, p2: Path):
    return Path(os.path.relpath(str(p1), str(p2)))

def make_subset_comparisons(inputs, output_dir, matches, sample_name_regex):
    for df_path in inputs:
        df = pd.read_pickle(df_path)
        if all(
                matches(
                    sample_name_regex.search(
                        Path(df[x + "sample"].iloc[0]).name
                    ).group(1)                    
                )
                for x in ["q", "s"]
        ):
            dest = output_dir / df_path.name
            dest.symlink_to(relative_to(df_path, dest.parent))
            yield df

def handle_filters(include, include_file):
    include = set(include)
    try:
        with open(include_file, "r") as filter_file:
            include |= set(l.rstrip() for l in filter_file)
    except TypeError:
        pass
    return include

def main():
    args = handle_arguments()
    args.output_dir.mkdir(exist_ok=True)
    matches = matcher(
        handle_filters(args.filter, args.filter_file),
        args.filter_regex
    )
    consume(
        make_subset_comparisons(
            tqdm(args.inputs),
            args.output_dir,
            matches,
            args.sample_name_regex
        )
    )
            
if __name__ == "__main__":
    main()
