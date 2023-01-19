import argparse
import sys
import math
import re
import functools
from collections import defaultdict, Counter
from tqdm import tqdm
from blasting import BlastnSearch
from heapq import nlargest

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
    parser.add_argument("--gene-group", "-g", type=int, default=1)
    parser.add_argument("-e", "--evalue", type=float, default=1e-50)
    parser.add_argument(
        "-n",
        "--top-n",
        type=int,
        default=1,
        help="use top n matches"
    )
    return parser.parse_args()

def parse_seq_id(regex, group, seq_id):
    #print(seq_id)
    return regex.match(seq_id).group(group)

def top(l, n=1, key=lambda x: x, value=lambda x: x):
    if n == -1:
        return list(l)
    value_counts = Counter(key(v) for v in l)
    top_keys = set(nlargest(n, value_counts))
    return [value(x) for x in l if key(x) in top_keys]
    # top_values = []
    # top_key = -math.inf
    # for v in l:
    #     k = key(v)
    #     if k > top_key:
    #         top_values = [value(v)]
    #         top_key = k
    #     elif k == top_key:
    #         top_values.append(value(v))
    # return top_values

def ind(x):
    def ind_inner(t):
        return t[x]
    return ind_inner

def gene_matches(parse, path1, path2, evalue, n=1):
    matches = defaultdict(list)
    for hit in tqdm(BlastnSearch(path1, path2, evalue=evalue)):
        matches[
            parse(hit.query_acc)
        ].append((hit.bit_score, parse(hit.subject_acc)))
    return {
        (k, su): sc
        for (k, l) in matches.items()
        for (sc, su) in top(l, key=ind(0), n=n)
    }

def tuprev(t):
    return tuple(reversed(t))

eprint = functools.partial(print, file=sys.stderr)

def main():
    args = parse_arguments()
    parse = functools.partial(parse_seq_id, args.regex, args.gene_group)
    gm = functools.partial(
        gene_matches,
        parse=parse,
        evalue=args.evalue,
        n=args.top_n
    )
    eprint("Getting forward matches.")
    forward_matches = gm(path1=args.transcripts1, path2=args.transcripts2)
    eprint("Getting backward matches.")
    backward_matches = {
        (b, a): c
        for ((a, b), c)
        in gm(path1=args.transcripts2, path2=args.transcripts1).items()
    }
    matches = set(forward_matches).intersection(set(backward_matches))
    match_lists = defaultdict(list)
    for (a, b) in matches:
        match_lists[a].append(
            (
                max(
                    forward_matches[(a, b)],
                    backward_matches[(a, b)]
                ),
                b
            )
        )
    matches = set(
        (a, b)
        for (a, l) in match_lists.items()
        for b in top(l, key=ind(0), value=ind(1))
    )
    for match in matches:
        print(*match)
    eprint(f"Found {len(matches)} matches.")
    
if __name__ == "__main__":
    main()
