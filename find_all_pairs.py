import argparse
import re
import multiprocessing
import functools
import itertools
from pathlib import Path
from find_homologs import HomologFinder, eprint
from blastdb_cache import BlastDBCache
from joblib import Parallel, delayed
from tqdm import tqdm

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputs", nargs="+", type=Path, required=True)
    parser.add_argument("-O", "--output-dir", type=Path, required=True)
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
        default=re.compile("^(.*?)_.*$")
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
    return parser.parse_args()

def find_homologs_and_save(
        transcripts1,
        transcripts2,
        out_path,
        hf_args=None,
        hf_kwargs=None
):
    if hf_args is None:
        hf_args = []
    if hf_kwargs is None:
        hf_kwargs = {}
    finder = HomologFinder(*hf_args, **hf_kwargs)
    table = finder.get_match_table(transcripts1, transcripts2)
    table["ssample"] = str(transcripts1)
    table["qsample"] = str(transcripts2)
    table.to_pickle(out_path)
    return True

def make_output_path(dir_, t1, t2, regex=None, extension="pkl"):
    t1 = t1.stem
    t2 = t2.stem
    if regex:
        t1 = regex.match(t1).group(1)
        t2 = regex.match(t2).group(1)
    return dir_ / ("{}--{}.{}".format(t1, t2, extension))

def make_one_db(db_loc, seq_file_path):
    cache = BlastDBCache(db_loc)
    cache.makedb(seq_file_path)
    return cache

def make_all_dbs(db_loc, seqs, jobs=1):
    cdict = {}
    for cache in Parallel(n_jobs=jobs)(
            delayed(make_one_db)(db_loc, p) for p in tqdm(seqs)
    ):
        cdict |= cache._cache
    cache = BlastDBCache(db_loc)
    cache._cache = cdict
    return cache

def main():
    args = handle_arguments()
    args.output_dir.mkdir(exist_ok=True)
    cache = None
    if args.db_cache_dir:
        args.db_cache_dir.mkdir(exist_ok=True)
        eprint("Building BLAST DBs.")
        cache = make_all_dbs(args.db_cache_dir, args.inputs, jobs=args.jobs)        
    fh = functools.partial(
        find_homologs_and_save,
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
    mop = functools.partial(
        make_output_path,
        regex=args.sample_regex,
        extension="pkl"
    )
    combos = list(itertools.combinations(args.inputs, 2))
    res = list(
        Parallel(n_jobs=args.jobs)(
            delayed(
                fh
            )(*p, mop(args.output_dir, *p))
            for p in tqdm(combos)
        )
    )
    assert all(res)
    
if __name__ == "__main__":
    main()
