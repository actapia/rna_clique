import argparse
import functools
import itertools

from git_archive_all import GitArchiver
from pathlib import Path

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("--extra-data", type=Path, nargs="*", default=[])
    return parser.parse_args()

def main():
    args = handle_arguments()
    GitArchiver(
        prefix=args.archive.stem,
        extra=itertools.chain.from_iterable(
            x.rglob("*") if x.is_dir() else x
            for x in args.extra_data
        )
    ).create(args.archive)
            
if __name__ == "__main__":
    main()

                
