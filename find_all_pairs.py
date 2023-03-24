import argparse
import re
import multiprocessing
import functools
import itertools
from pathlib import Path
from find_homologs import HomologFinder
from joblib import Parallel, delayed
from tqdm import tqdm

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputs", nargs="+", type=Path, required=True)
    parser.add_argument("-O", "--output-dir", type=Path, required=True)
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

def main():
    args = handle_arguments()
    args.output_dir.mkdir(exist_ok=True)
    fh = functools.partial(
        find_homologs_and_save,
        hf_args=[
            args.gene_regex,
            args.top_n,
            args.evalue,
            args.keep_all
        ]
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
