import sys
import csv
import argparse
import re
import sys
import os

from pathlib import Path

from more_itertools import windowed, padded

arg_re = re.compile("--?(.*)")
cache_threshold = 5

def handle_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ignore-positional", action="store_true")
    parser.add_argument(
        "-c",
        "--cache",
        nargs="?",
        choices=[True, False, "auto"],
        default=False,
        const=True,
    )
    return parser.parse_args()

def main():
    args = handle_args()
    if "OUTDIR_ROOT" in os.environ:
        root = Path(os.environ["OUTDIR_ROOT"])
    else:
        root = Path(".")
    writer = csv.writer(sys.stdout, delimiter=" ")
    pcount = 0
    for line in csv.reader(sys.stdin, delimiter=" "):
        kwargs = {}
        matches = [arg_re.match(v) for v in line]
        wind = windowed(
            padded(zip(line, matches), fillvalue=(None, None), n=len(line)+1),
            2
        )
        for (arg, match), (next_arg, next_match) in wind:
            if match:
                if not next_match and next_arg:
                    kwargs[match.group(1)] = next_arg
                    # skip = True
                    next(wind)
                else:
                    kwargs[match.group(1)] = ""
            else:
                pcount += 1
                if not args.ignore_positional:
                    kwargs[arg] = ""
        outdir = root / Path(
            "out_{}".format("_".join(k+v for (k, v) in kwargs.items()))
        )
        line = line + [
            "--out-dir-1", str(outdir / "od1"),
            "--out-dir-2", str(outdir / "od2"),
            "--output-graph", str(outdir / "graph.pkl")
        ]
        if args.cache is True or \
           (args.cache == "auto" and pcount > cache_threshold):
            line = line + ["--cache-dir", str(outdir / "dbcache")]
        writer.writerow(line)
        
if __name__ == "__main__":
    main()

