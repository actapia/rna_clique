import argparse
import sys

from pathlib import Path

import Bio.SeqIO
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("transcripts", type=Path)
    parser.add_argument("quant", type=Path)
    args = parser.parse_args()
    quant = pd.read_table(args.quant).set_index("Name")
    #from IPython import embed; embed()
    for transcript in Bio.SeqIO.parse(args.transcripts, "fasta"):
        transcript.description = "tpm{}_{}".format(
            quant.TPM[transcript.id],
            transcript.description
        )
        transcript.id = "tpm{}_{}".format(
            quant.TPM[transcript.id],
            transcript.id
        )
        transcript.name = ""
        Bio.SeqIO.write(transcript, sys.stdout, "fasta")

if __name__ == "__main__":
    main()
