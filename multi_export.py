import argparse
import sys
import os
from pathlib import Path
from export_orthologs import OrthologExporter
from filtered_distance import SampleSimilarity
from make_subset import multi_glob

from tqdm import tqdm

def common_suffix(words):
    # From yairchu on Stack Overflow
    # https://stackoverflow.com/a/76563105/1927609
    return os.path.commonprefix([w[::-1] for w in words])[::-1]

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-A",
        "--analyses",
        type=Path,
        nargs="+"
    )
    parser.add_argument(
        "-O",
        "--out-dir",
        type=Path,
        default=Path(".")
    )
    # parser.add_argument(
    #     "--zip",
    #     action="store_true"
    # )
    parser.add_argument(
        "--strip-common-suffix",
        "-s",
        action="store_true"
    )
    return parser.parse_args()

def main():
    args = handle_arguments()
    args.out_dir.mkdir(exist_ok=True)
    for components in range(1, max(len(a.parts) for a in args.analyses)):
        existing = set()
        for analysis in args.analyses:
            name = "_".join(analysis.parts[-components:])
            if name in existing:
                break
            existing.add(name)
        else:
            suffix = common_suffix(existing)
            break
    else:
        print(
            ("Could not automatically assign names to analysis outputs"
             "directories."),
            file=sys.stderr
        )
    # print(existing)
    # sys.exit(1)
    for analysis in tqdm(args.analyses):
        dir_name = "_".join(analysis.parts[-components:])
        if args.strip_common_suffix:
            dir_name = dir_name.removesuffix(suffix)
        out_dir = args.out_dir / dir_name
        out_dir.mkdir(exist_ok=True)
        comparisons = list(multi_glob(analysis / "od2", ["*.pkl", "*.h5"]))
        graph = analysis / "graph.pkl"
        sim = SampleSimilarity.from_filenames(
            graph,
            comparisons,
            store_dfs=True
        )

        for i, name in enumerate(["contributing", "all"]):
            sub_dir = out_dir / name
            sub_dir.mkdir(exist_ok=True)
            exporter = OrthologExporter(
                sim,
                non_contributing=i,
            )
            exporter.by_component(sub_dir, order="before")

            
    
if __name__ == "__main__":
    main()
