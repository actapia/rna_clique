import os
import ordered_sample_count_test
from genotype_endophyte_interleave_test import handle_arguments
from load_sample_metadata import get_sample_metadata
from fair_sample_count_tests import make_arg_arr

from IPython import embed


def main():
    args, remainder = handle_arguments()
    extra_args = make_arg_arr(remainder)
    sample_order = get_sample_metadata(args.sample_metadata)["name"]
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
    main()
