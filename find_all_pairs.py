import argparse
import re
import math
import multiprocessing
import functools
import itertools
import config as config_module
from more_itertools import consume
from typing import Optional, Any, Callable
from collections.abc import Iterable, Mapping
from pathlib import Path
from find_homologs import HomologFinder, eprint
from simple_blast import BlastDBCache
from joblib import Parallel, delayed
import os
from tqdm import tqdm

from gene_matches_tables import write_table
from transcripts import TranscriptID

default_sample_regex = re.compile(os.environ.get("SAMPLE_RE", "^(.*?)_.*$"))

def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager()
    arg_config.expose_fields_with_default_aliases(
        "top_genes_dir",
        "output_dir",
        "transcript_id_regex",
        required=False
    )
    arg_config.expose_fields_with_default_aliases(
        "cache_dir",
        "evalue",
        "title",
    )
    arg_config.add_argument(
        "--sample-regex",
        "-R",
        type=re.compile,
        default=default_sample_regex,
        help="Python regex for parsing sample names"
    )
    arg_config.add_output_config_argument()
    return arg_config

def find_homologs_and_save(
        transcripts1 : Path,
        transcripts2 : Path,
        out_path : Path,
        hf_args : Optional[Iterable] = None,
        hf_kwargs : Optional[Mapping[str, Any]] = None
):
    """Get the gene matches tables for the given FASTA files and save results.

    Parameters:
        transcripts1: Path to the top n transcripts FASTA for the first sample.
        transcripts2: Path to the top n transcripts FASTA for the second sample.
        out_path:     Output file in which to store the gene matches table.
        hf_args:      Arguments to pass to HomologFinder constructor.
        hf_kwargs:    Keyword arguments to pass to HomologFinder constructor.
    """
    if hf_args is None:
        hf_args = []
    if hf_kwargs is None:
        hf_kwargs = {}
    finder = HomologFinder(*hf_args, **hf_kwargs)
    table = finder.get_match_table(transcripts1, transcripts2)
    table["ssample"] = str(transcripts1)
    table["qsample"] = str(transcripts2)
    table[["ssample", "qsample"]] = table[["ssample", "qsample"]].astype(
        "category"
    )
    write_table(table, out_path)
    return table

def make_output_path(
        dir_ : Path,
        t1 : Path,
        t2 : Path,
        path_to_sample: Optional[Callable] = None,
        extension : str = "h5"
) -> Path:
    """Return the output path for the comparison between the two files.

    Parameters:
        dir_:            Path to output directory.
        t1:              Path to (top n) transcripts for first sample.
        t2:              Path to (top n) transcripts for second sample.
        path_to_sample:  Function mapping paths to sample names.
        extension (str): File extension to use for output file.

    Returns:
        The constructed path for the comparison between t1 and t2.
    """
    # t1 = t1.stem
    # t2 = t2.stem
    if path_to_sample:
        t1 = path_to_sample(t1)
        t2 = path_to_sample(t2)
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
        seqs:       Paths to the FASTA files for which to make the databases.
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

def find_all_pairs(
        inputs,
        output_dir,
        cache_dir,
        path_to_sample,
        hf_args=[],
        jobs=multiprocessing.cpu_count()-1,
):
    cache = None
    if cache_dir:
        eprint("Building BLAST DBs.")
        cache = make_all_dbs(cache_dir, inputs, jobs=jobs)
    fh = functools.partial(
        find_homologs_and_save,
        hf_args=hf_args,
        hf_kwargs = {
            "db_cache": cache
        }
    )
    mop = functools.partial(
        make_output_path,
        path_to_sample=path_to_sample,
        extension="h5"
    )
    return (
        Parallel(n_jobs=jobs, return_as="generator_unordered")(
            delayed(
                fh
            )(*p, mop(output_dir, *p))
            for p in itertools.combinations(inputs, 2)
        ), map(
            lambda x: mop(output_dir, *x),
            itertools.combinations(inputs,2)
        ), math.comb(len(inputs), 2)
    )

def main():
    original_args, args, config = build_parser().get_arguments_and_config()
    config_module.RNACliqueConfigArgumentManager.make_output_dirs(config)
    id_parser = TranscriptID.parser_from_re(config.transcript_id_regex)
    if original_args.sample_regex or not config.path_to_sample:
        path_to_sample = lambda x: args.sample_regex.match(x.stem).group(1)
    else:
        path_to_sample = config.path_to_sample            
    gen, _, gen_len = find_all_pairs(
        args.inputs,
        args.output_dir,
        args.db_cache_dir,
        path_to_sample,
        hf_args=[
            id_parser,
            args.top_n,
            args.evalue,
            args.keep_all
        ]
    )
    consume(tqdm(gen, total=gen_len))
    config.mark_finish()
    config.yaml_save(args.output_config)
    
if __name__ == "__main__":
    main()
