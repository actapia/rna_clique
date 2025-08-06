import heapq
import argparse
import re
import sys

import Bio
import Bio.SeqIO

from pathlib import Path
from collections import defaultdict
from collections.abc import Collection

from typing import Iterator, Callable

default_gene_re = re.compile(r"^.*cov_([0-9]+(?:\.[0-9]+))_g([0-9]+)_i([0-9]+)")

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--transcripts", "-i", type=Path)
    parser.add_argument("--top", "-n", type=int, required=True)
    parser.add_argument(
        "--pattern",
        "-e",
        type=re.compile,
        default=default_gene_re
    )
    return parser.parse_args()

class TopGeneSelector:
    def __init__(
            self,
            transcripts: Callable[[], Iterator[Bio.SeqRecord]],
            top: int,
            gene_re: re.Pattern
    ):
        self.transcripts = transcripts
        self.top = top
        self.gene_re = gene_re

    def get_top_genes(self):
        highest_coverage = defaultdict(float)
        for t in self.transcripts():
            cov, gene, iso = self.gene_re.match(t.id).groups()
            gene = int(gene)
            highest_coverage[gene] = max(highest_coverage[gene], float(cov))
        for _, k in heapq.nlargest(
                self.top,
                ((v, k) for (k, v) in highest_coverage.items())
        ):
            yield k

    def get_top_gene_seqs(self):
        top_genes = set(self.get_top_genes())
        for t in self.transcripts():
            cov, gene, iso = self.gene_re.match(t.id).groups()
            if int(gene) in top_genes:
                yield t

    @classmethod
    def from_path(cls, path: Path, *args, **kwargs):
        return cls(lambda: Bio.SeqIO.parse(path, "fasta"), *args, **kwargs)

    @classmethod
    def from_sequences(cls, seqs: Collection[Bio.SeqRecord], *args, **kwargs):
        return cls(lambda: seqs, *args, **kwargs)
        
def main():
    args = handle_arguments()
    if args.transcripts:
        top = TopGeneSelector.from_path(
            args.transcripts,
            args.top,
            args.pattern
        )
    else:
        top = TopGeneSelector.from_sequences(
            list(Bio.SeqIO.parse(sys.stdin, "fasta")),
            args.top,
            args.pattern
        )
    Bio.SeqIO.write(top.get_top_gene_seqs(), sys.stdout, "fasta")

if __name__ == "__main__":
    main()

