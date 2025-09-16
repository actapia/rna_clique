import heapq
import sys

import Bio
import Bio.SeqIO

from pathlib import Path
from collections import defaultdict
from collections.abc import Collection
from typing import Iterator, Callable

from .transcripts import TranscriptID
from . import config as config_module

def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager()
    arg_config.expose_fields_with_default_aliases(
        "top_genes",
        "transcript_id_regex",
        required=True
    )
    arg_config.add_argument("--transcripts", "-i", type=Path)
    return arg_config

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
    _, args, config = build_parser().get_arguments_and_config()
    parse_transcript_id = TranscriptID.parser_from_re(
        config.transcript_id_regex
    )
    if args.transcripts:
        top = TopGeneSelector.from_path(
            args.transcripts,
            config.top_genes,
            parse_transcript_id
        )
    else:
        top = TopGeneSelector.from_sequences(
            list(Bio.SeqIO.parse(sys.stdin, "fasta")),
            config.top_genes,
            parse_transcript_id,
        )
    Bio.SeqIO.write(top.get_top_gene_seqs(), sys.stdout, "fasta")

if __name__ == "__main__":
    main()

