import argparse
import re
import multiprocessing
import functools
import itertools
from typing import Optional
from collections.abc import Iterable
from pathlib import Path
from find_homologs import HomologFinder, eprint
from simple_blast import BlastDBCache
from joblib import Parallel, delayed
import os
from tqdm import tqdm

default_sample_regex = re.compile(os.environ.get("SAMPLE_RE", "^(.*?)_.*$"))

def handle_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "get comparisons (gene matches tables) for all pairs of "
            "FASTA files"
        )
    )
    parser.add_argument(
        "-i",
        "--inputs",
        nargs="+",
        type=Path,
        required=True,
        help="input FASTA files (containing top n genes)"
    )
    parser.add_argument(
        "-O",
        "--output-dir",
        type=Path,
        required=True,
        help="directory in which to write BLAST results"
    )
    parser.add_argument(
        "-D",
        "--db-cache-dir",
        type=Path,
        help="directory in which to store BLAST DBs for input FASTA files"
    )
    parser.add_argument(
        "--gene-regex",
        "-r",
        type=re.compile,
        default=re.compile("^.*g([0-9]+)_i([0-9]+)"),
        help="Python regex for parsing sequence IDs"
    )
    parser.add_argument(
        "--sample-regex",
        "-R",
        type=re.compile,
        default=default_sample_regex,
        help="Python regex for parsing sample names"
    )
    parser.add_argument(
        "-e",
        "--evalue",
        type=float,
        default=1e-50,
        help="e-value threshold to use for BLAST alignemnts"
    )
    parser.add_argument(
        "-n",
        "--top-n",
        type=int,
        default=1,
        help="use top n matches (parameter big N)"
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
        default=multiprocessing.cpu_count()-1,
        help="number of parallel jobs to use"
    )
    return parser.parse_args()

def find_homologs_and_save(
        transcripts1 : Path,
        transcripts2 : Path,
        out_path : Path,
        hf_args=None,
        hf_kwargs=None
):
    """Get the gene matches tables for the given FASTA files and save results.

    Parameters:
        transcripts1: Path to the top n transcripts FASTA for the first sample.
        transcripts2: Path to the top n transcripts FASTA for the second sample.
        out_path:     Output file in which to store the gene matches table.
    """
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

def make_output_path(
        dir_ : Path,
        t1 : Path,
        t2 : Path,
        regex : Optional[re.Pattern] = None,
        extension : str = "pkl"
) -> Path:
    """Return the output path for the comparison between the two files.

    Parameters:
        dir_:            Path to output directory.
        t1:              Path to (top n) transcripts for first sample.
        t2:              Path to (top n) transcripts for second sample.
        regex:           Regular expression to use to parse t1 and t2 filenames.
        extension (str): File extension to use for output file.

    Returns:
        The constructed path for the comparison between t1 and t2.
    """
    t1 = t1.stem
    t2 = t2.stem
    if regex:
        t1 = regex.match(t1).group(1)
        t2 = regex.match(t2).group(1)
    return dir_ / ("{}--{}.{}".format(t1, t2, extension))

def make_one_db(db_loc : Path, seq_file_path : Path) -> BlastDBCache:
    """Create a BlastDBCache with a database for a single FASTA file.

    Parameters:
        db_loc:        The directory in which to create the BLAST DB.
        seq_file_path: The path to the subject sequence FASTA.

    Returns:
        A BlastDBCache with a DB for the provided FASTA file at seq_file_path.
    """
    cache = BlastDBCache(db_loc)
    cache.makedb(seq_file_path)
    return cache

def make_all_dbs(
        db_loc : Path,
        seqs : Iterable[Path],
        jobs : int = 1
) -> BlastDBCache:
    """Create a BlastDBCache with databases for the given FASTA files.

    Since one database does not depend on any other, this function optionally
    creates the databases in parallel.

    Parameters:
        db_loc:     Path to the directory in which to make the databases.
        seqs:       Paths to the FASTA files for which to makes the databases.
        jobs (int): The number of parallel jobs to use.

    Returns:
        A BlastDBCache with databases for all provided FASTA files.
    """
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
