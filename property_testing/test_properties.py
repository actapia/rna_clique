import argparse
import sys
import re
import itertools
import functools
import multiprocessing
from fractions import Fraction
from collections import defaultdict
from pathlib import Path

from joblib import Parallel, delayed
from tqdm import tqdm

from find_all_pairs import make_all_dbs, default_sample_regex

from find_homologs import HomologFinder, eprint

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputs", nargs="+", type=Path, required=True)
    parser.add_argument("-D", "--db-cache-dir", type=Path)
    parser.add_argument(
        "--gene-regex",
        "-r",
        type=re.compile,
        default=re.compile("^.*g([0-9]+)_i([0-9]+)")
    )
    parser.add_argument(
        "--sample-regex",
        "-R",
        type=re.compile,
        default=default_sample_regex
    )
    parser.add_argument("-e", "--evalue", type=float, default=1e-50)
    parser.add_argument(
        "-n",
        "--top-n",
        type=int,
        default=1,
        help="use top n matches"
    )
    parser.add_argument(
        "--keep-all",
        "-k",
        action="store_true",
        help="keep all pairs in case of a tie"
    )
    parser.add_argument(
        "--jobs",
        "-j",
        type=int,
        default=multiprocessing.cpu_count()-1
    )
    parser.add_argument(
        "--debug",
        action="store_true"
    )
    return parser.parse_args()

def large_sum(col):
    try:
        return sum(r.item() for r in col)
    except AttributeError:
        return sum(r for r in col)

def find_homologs(
        transcripts1,
        transcripts2,
        hf_args=None,
        hf_kwargs=None
):
    if hf_args is None:
        hf_args = []
    if hf_kwargs is None:
        hf_kwargs = {}
    finder = HomologFinder(*hf_args, **hf_kwargs)
    table = finder.get_match_table(transcripts1, transcripts2)
    return transcripts1, transcripts2, Fraction(
        large_sum(table["nident"]),
        large_sum(table["length"]) - large_sum(table["gaps"])
    )

def main():
    args = handle_arguments()
    cache = None
    if args.db_cache_dir:
        args.db_cache_dir.mkdir(exist_ok=True)
        eprint("Building BLAST DBs.")
        cache = make_all_dbs(args.db_cache_dir, args.inputs, jobs=args.jobs)
    perms = list(itertools.permutations(args.inputs, 2))
    distances = defaultdict(dict)
    fh = functools.partial(
        find_homologs,
        hf_args=[
            args.gene_regex,
            args.top_n,
            args.evalue,
            args.keep_all
        ],
        hf_kwargs={
            "db_cache": cache
        }
    )
    fail = False
    for t1, t2, dist in Parallel(n_jobs=args.jobs, return_generator=True)(
            delayed(
                fh
            )(*p)
            for p in tqdm(perms)
    ):
        if distances[t2].get(t1, dist) != dist:
            eprint(
                "{} depends on order ({} != {})".format(
                    {t1, t2},
                    distances[t2][t1],
                    dist
                )
            )
            fail = True
        elif args.debug and t1 in distances[t2]: 
            eprint(
                "\033[92m{} does not depend on order ({} == {})\033[0m".format(
                    {t1, t2},
                    distances[t2][t1],
                    dist
                )
            )
        distances[t1][t2] = dist
    # Check triangle inequality.
    for t1, t2, t3 in itertools.permutations(args.inputs, 3):
        d1 = 1 - distances[t1][t2]
        d2 = 1 - distances[t1][t3]
        d3 = 1 - distances[t3][t2]
        if d1 > d2 + d3:
            eprint(
                ("{} does not satisfy the triangle inequality "
                 "({} > {} + {})").format(
                     [str(t) for t in (t1, t2, t3)],
                     d1,
                     d2,
                     d3
                 )
            )
            fail = True
    return fail

if __name__ == "__main__":
    sys.exit(main())



# from find_homologs import HomologFinder

# def handle_arguments():
#     parser = argparse.ArgumentParser()
#     parser.add_argument("--transcripts", "-t", nargs="+", type=Path)
#     parser.add_argument(
#         "--regex",
#         "-r",
#         type=re.compile,
#         default=re.compile("^.*g([0-9]+)_i([0-9]+)")
#     )
#     parser.add_argument("-e", "--evalue", type=float, default=1e-50)
#     parser.add_argument(
#         "-n",
#         "--top-n",
#         type=int,
#         default=1,
#         help="use top n matches"
#     )
#     args = parser.parse_args()
#     return args.parse_args()

# def main():
#     args = handle_arguments()
#     for s1, s2 in combinations(args.transcripts, 2):



# if __name__ == "__main__":
#     main()
