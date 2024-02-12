import argparse
#import re
import os
import ordered_sample_count_test
import pandas as pd
from pathlib import Path
from load_sample_metadata import get_sample_metadata
from fair_sample_count_tests import make_arg_arr

from IPython import embed

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-metadata", "-s", type=Path, required=True)
    return ordered_sample_count_test.handle_arguments(parser)

def group_interleave(df, cols):
    return pd.concat(
        t[1].reset_index(drop=True)
        for t in df.groupby(cols)
    ).sort_index().reset_index(drop=True)

def do_group_interleave_test(cols):
    args, remainder = handle_arguments()
    extra_args = make_arg_arr(remainder)
    sample_info = get_sample_metadata(args.sample_metadata)
    sample_order = group_interleave(sample_info, cols)["name"]
    fasta_dict = {
        os.path.basename(f): f
        for f in args.fasta
    }
    ordered_sample_count_test.OrderedSampleCountTest(
        (fasta_dict[f] for f in sample_order),
        args.out_dir,
        args=extra_args
    ).run()

if __name__ == "__main__":
    do_group_interleave_test(["genotype", "endophyte"])
    

    
