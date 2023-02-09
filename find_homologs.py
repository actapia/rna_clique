import argparse
import sys
import re
import functools
from blasting import BlastnSearch
import numpy as np

import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("transcripts1")
    parser.add_argument("transcripts2")
    parser.add_argument(
        "--regex",
        "-r",
        type=re.compile,
        default=re.compile("^.*g([0-9]+)_i([0-9]+)")
    )
    parser.add_argument("-e", "--evalue", type=float, default=1e-50)
    parser.add_argument(
        "-n",
        "--top-n",
        type=int,
        default=1,
        help="use top n matches"
    )
    return parser.parse_args()

def parse_seq_id(regex, s):
    return s.str.extract(regex).astype(np.int32)

def gene_matches(parse, path1, path2, evalue, n=1):
    search = BlastnSearch(path1, path2, evalue=evalue)
    search.hits[["qgene","qiso"]] = parse(search.hits["qseqid"])
    search.hits[["sgene","siso"]] = parse(search.hits["sseqid"])
    return highest_bitscores(search.hits, n, keep="all")

eprint = functools.partial(print, file=sys.stderr)

def highest_bitscores(df, n=1, groupby="qgene", **kwargs):
    return df.loc[
        df.groupby(
            groupby
        )["bitscore"].nlargest(n, **kwargs).index.get_level_values(-1)
    ]

def main():
    args = parse_arguments()
    parse = functools.partial(parse_seq_id, args.regex)
    gm = functools.partial(
        gene_matches,
        parse=parse,
        evalue=args.evalue,
        n=args.top_n
    )
    eprint("Getting forward matches.")
    forward_matches = gm(path1=args.transcripts1, path2=args.transcripts2)
    eprint("Getting backward matches.")
    backward_matches = gm(path1=args.transcripts2, path2=args.transcripts1)
    backward_matches.rename(
        columns={
            a+v : b+v
            for (a,b) in [("q","s"),("s","q")]
            for v in ["seqid","gene","iso"]
        },
        inplace=True
    )
    merge_columns = ["qgene", "sgene"]
    intersection = pd.merge(
        forward_matches[merge_columns].reset_index(),
        backward_matches[merge_columns].reset_index(),
        how="inner",
        on=merge_columns
    )
    best_matches = highest_bitscores(
        highest_bitscores(
            pd.concat(
                [
                    forward_matches.loc[
                        intersection["index_x"].drop_duplicates()
                    ],
                    backward_matches.loc[
                        intersection["index_y"].drop_duplicates()
                    ],
                ]
            ).reset_index(drop=True),
            groupby=merge_columns
        ),
    )    
    dedup = best_matches[merge_columns].drop_duplicates()
    for match in dedup.itertuples(index=False):
        print(*match)
    eprint(f"Found {len(dedup)} matches.")
    
if __name__ == "__main__":
    main()
