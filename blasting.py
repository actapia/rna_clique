import subprocess
import pandas as pd

default_out_columns = ['qseqid',
 'sseqid',
 'pident',
 'length',
 'mismatch',
 'gapopen',
 'qstart',
 'qend',
 'sstart',
 'send',
 'evalue',
 'bitscore']

class BlastnSearch:
    def __init__(
            self,
            seq1_path,
            seq2_path,
            evalue=1e-20,
            out_columns=default_out_columns,
            additional_columns=[]
    ):
        self._seq1_path = seq1_path
        self._seq2_path = seq2_path
        self._evalue = evalue
        self._hits = None
        self._out_columns = list(out_columns + additional_columns)

    @property
    def seq1_path(self):
        return self._seq1_path

    @property
    def seq2_path(self):
        return self._seq2_path

    @property
    def evalue(self):
        return self._evalue

    @property
    def hits(self):
        if self._hits is None:
            self._get_hits()
        return self._hits

    def __len__(self):
        return len(self.hits)

    def __iter__(self):
        yield from self.hits

    def _get_hits(self):
        proc = subprocess.Popen(
            [
                "blastn",
                "-subject",
                self.seq1_path,
                "-query",
                self.seq2_path,
                "-evalue",
                str(self.evalue),
                "-outfmt",
                " ".join(["6"] + self._out_columns)
            ]
            ,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        #self._hits = tuple(BlastHit.parse_blast_lines(out.decode("utf-8").splitlines()))
        # for l in proc.stderr:
        #     print(l)
        self._hits = pd.read_csv(proc.stdout, names=self._out_columns, delim_whitespace=True)
