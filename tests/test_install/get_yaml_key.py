from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

import argparse
from pathlib import Path

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=Path)
    parser.add_argument("key", type=str)
    return parser.parse_args()

def main():
    args = handle_arguments()
    with open(args.path, "r") as f:
        print(load(f, Loader=Loader)[args.key])

if __name__ == "__main__":
    main()
