import argparse
import pickle
import re
from filtered_distance import SampleSimilarity
from path_to_sample import path_to_sample

from pathlib import Path

import Bio.SeqIO

from tqdm import tqdm
#from IPython import embed

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--graph", type=Path, required=True)
    parser.add_argument(
        "-c",
        "--comparisons",
        type=Path,
        nargs="+",
        required=True
    )
    parser.add_argument("--rename", action="store_true")
    parser.add_argument(
        "-s",
        "--samples",
        type=int
    )
    parser.add_argument(
        "--gene-regex",
        "-r",
        type=re.compile,
        default=re.compile("^.*g([0-9]+)_i([0-9]+)"),
        help="Python regex for parsing sequence IDs"
    )
    parser.add_argument(
        "--out-dir",
        "-O",
        type=Path,
        required=True
    )
    return parser.parse_args()

def main():
    args = handle_arguments()
    args.out_dir.mkdir(exist_ok=True)
    sim = SampleSimilarity.from_filenames(
        args.graph,
        tqdm(args.comparisons),
        sample_count=args.samples
    )
    samples = sim.valid["sample"].drop_duplicates()
    valid_tuples = set(map(tuple, sim.valid.itertuples(index=False)))
    for sample in tqdm(samples):
        out_fn = args.out_dir / "{}_orthologs.fasta".format(
            path_to_sample(sample)
        )
        gen = (
                seq for seq in Bio.SeqIO.parse(sample, "fasta")
                if (
                            sample,
                            int(args.gene_regex.search(seq.id).group(1))
                ) in valid_tuples
            )
        Bio.SeqIO.write(
            gen,
            out_fn,
            "fasta"
        )

    
if __name__ == "__main__":
    main()
