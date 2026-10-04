import sys
import re

from collections import defaultdict

import Bio.SeqIO

def counter():
    """Returns a function returning and incrementing counter i on each call."""
    def increment():
        nonlocal i
        i += 1
        return i
    i = -1
    return increment

trinity_gene_id_re = re.compile(r"TRINITY_DN\d+_c\d+_g\d+")

def main():
    ids = defaultdict(counter())
    for transcript in Bio.SeqIO.parse(sys.stdin, "fasta"):
        gene_match = trinity_gene_id_re.search(transcript.id)
        new_id = "{}_gid{}".format(
            gene_match.group(0),
            ids[gene_match.group(0)]
        )
        transcript.id = trinity_gene_id_re.sub(new_id, transcript.id)
        transcript.description = trinity_gene_id_re.sub(
            new_id,
            transcript.description
        )
        Bio.SeqIO.write(transcript, sys.stdout, "fasta")

if __name__ == "__main__":
    main()
