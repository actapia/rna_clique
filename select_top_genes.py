import heapq
import sys

import Bio
import Bio.SeqIO

import config as config_module

from pathlib import Path
from collections import defaultdict
from collections.abc import Collection
from typing import Iterator, Callable

from transcripts import TranscriptID, default_parser

def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager(
        description=(
            "Select top genes by k-mer coverage in a transcripts FASTA file."
        )
    )
    arg_config.expose_fields_with_default_aliases(
        "top_genes",
        "transcript_id_regex",
        required=True
    )
    arg_config.add_argument(
        "--transcripts",
        "-i",
        type=Path,
        help="FASTA file from which to select top genes by k-mer coverage.",
    )
    return arg_config

class TopGeneSelector:
    """Class for selecting top n genes by k-mer coverage from transcripts.

    This class contains methods for obtaining the top n genes in a transcriptome
    according to their k-mer coverage. k-mer coverage is used here to roughly
    measure the extent to which genes were expressed in the RNA-seq libraries
    from which the transcriptomes were derived.

    Since k-mer coverages are defined for individual transcripts rather than
    genes, this class defines the k-mer coverage of a gene as the maximum k-mer
    coverage among the gene's isoforms.

    Note that this class always selects exactly n genes. It does NOT find the
    top n distinct gene k-mer coverages and then select the genes that have
    those coverages. Instead, it behaves as though it sorts the genes by k-mer
    coverage in descending order, and then selects the first n elements of that
    sorted collection. Among other things, this means that tie" may be broken
    arbitrarily; it is possible that the last element selected has the same
    k-mer coverage as another element that was not selected.

    Efficiently obtaining both the IDs of the top n genes and the sequences of
    all transcript isoforms associated with those genes requires iterating over
    the sequences in the provided transcriptome twice---once for each of those
    two steps. Since sometimes it may only be necessary to obtain the IDs, and
    the sequences themselves can be ignored, this class contains two
    methods---one to obtain just the gene IDs and one for the actual sequences.

    How the transcripts should best be iterated twice depends on the application
    and the available hardware resources. A typically slower, but less memory
    intensive, way of iterating twice is to simply parse the input FASTA file
    twice. In principle, this can be done for a FASTA file stored on disk
    without opening the file twice; fseek could be used instead. However,
    BioPython's Bio.SeqIO.parse function does not work well with seeking, so it
    is usually easier to simply call parse twice, opening the file each
    time. For applications and environments where memory is less of a concern,
    it may be faster to store all of the transcript sequences in memory (e.g.,
    in a list) and simply iterate over the in-memory collection twice. This
    class tries to accommodate both cases (and possibly others) by allowing the
    client code to pass a nullary Callable returning an iterator over
    transcripts. Since this makes the interface somewhat cumbersome for common
    use cases, the class provides classmethods that construct TopGeneSelector
    objects with such Callables automatically.

    Attributes:
        transcripts:         Function returning transcript SeqRecord iterator.
        top (int):           Number of top genes to select.
        parse_transcript_id: Function to parse FASTA IDs into TranscriptIDs.
    """
    def __init__(
            self,
            transcripts: Callable[[], Iterator[Bio.SeqRecord]],
            top: int,
            parse_transcript_id: Callable[[str], TranscriptID] = default_parser,
    ):
        """Construct a TopGeneSelector for given transcripts and top gene count.

        The first argument, transcripts, must be a nullary function returning an
        iterator over the transcript Bio.SeqRecrod objects from which to select
        the top genes. For convenience, this method also provides from_
        classmethods. Each can construct a TopGeneSelector object from a more
        common object, such as a Path to the FASTA file containing the
        transcripts or a list of Bio.SeqRecord objects.

        Parameters:
            transcripts:         Function to get transcript SeqRecord iterator.
            top (int):           Number of top genes to select.
            parse_transcript_id: Function to parse FASTA IDs into TranscriptIDs.
        """
        self.transcripts = transcripts
        self.top = top
        self.parse_transcript_id = parse_transcript_id

    def get_top_genes(self) -> Iterator[int]:
        """Get the gene IDs of the top genes by k-mer coverage."""
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
        """Get the Bio.SeqRecord objects of the top genes by k-mer coverage."""
        top_genes = set(self.get_top_genes())
        for t in self.transcripts():
            cov, gene, iso = self.parse_transcript_id(t.id)
            if int(gene) in top_genes:
                yield t

    @classmethod
    def from_path(cls, path: Path, *args, **kwargs):
        """Get a TopGeneSelector from a Path to the transcripts FASTA file."""
        return cls(lambda: Bio.SeqIO.parse(path, "fasta"), *args, **kwargs)

    @classmethod
    def from_sequences(cls, seqs: Collection[Bio.SeqRecord], *args, **kwargs):
        """Get a TopGeneSelector for a Collection of Bio.SeqRecord objects."""
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

