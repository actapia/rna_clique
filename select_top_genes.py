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

from transcripts import default_gene_re, TranscriptID

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
            parse_transcript_id: Callable[[str], TranscriptID]
    ):
        self.transcripts = transcripts
        self.top = top
        self.parse_transcript_id = parse_transcript_id

    def get_top_genes(self):
        highest_coverage = defaultdict(float)
        for t in self.transcripts():
            cov, gene, iso = self.parse_transcript_id(t.id)
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
            cov, gene, iso = self.parse_transcript_id(t.id)
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
    parse_transcript_id = TranscriptID.parser_from_re(args.pattern)
    if args.transcripts:
        top = TopGeneSelector.from_path(
            args.transcripts,
            args.top,
            parse_transcript_id
        )
    else:
        top = TopGeneSelector.from_sequences(
            list(Bio.SeqIO.parse(sys.stdin, "fasta")),
            args.top,
            parse_transcript_id,
        )
    Bio.SeqIO.write(top.get_top_gene_seqs(), sys.stdout, "fasta")

if __name__ == "__main__":
    main()

