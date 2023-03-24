import sys
import csv
import argparse
import re

from pathlib import Path

from more_itertools import windowed, padded

arg_re = re.compile("--?(.*)")

def handle_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ignore-positional", action="store_true")
    return parser.parse_args()

def main():
    args = handle_args()
    writer = csv.writer(sys.stdout, delimiter=" ")
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
            elif not args.ignore_positional:
                kwargs[arg] = ""
        outdir = Path(
            "out_{}".format("_".join(k+v for (k, v) in kwargs.items()))
        )
        line = line + [
            "--out-dir-1", str(outdir / "od1"),
            "--out-dir-2", str(outdir / "od2"),
            "--output-graph", str(outdir / "graph.pkl")
        ]
        writer.writerow(line)
        
if __name__ == "__main__":
    main()

