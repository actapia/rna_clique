import subprocess

class BlastnSearch:
    def __init__(self, seq1_path, seq2_path, evalue=1e-20):
        self._seq1_path = seq1_path
        self._seq2_path = seq2_path
        self._evalue = evalue
        self._hits = None

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
                "6"
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        out, err = proc.communicate()
        self._hits = tuple(BlastHit.parse_blast_lines(out.decode("utf-8").splitlines()))

# class BlastHits:
#     def __init__(self, hits):
#         self._hits = hits

#     @property
#     def 

class BlastHit:
    _field_types = [
        str,
        str,
        float,
        int,
        int,
        int,
        int,
        int,
        int,
        int,
        float,
        float
    ]
    
    def __init__(self,
                 query_acc,
                 subject_acc,
                 percent_identity,
                 alignment_length,
                 mismatches,
                 gap_opens,
                 query_start,
                 query_end,
                 subject_start,
                 subject_end,
                 evalue,
                 bit_score):
        self.query_acc = query_acc
        self.subject_acc = subject_acc
        self.percent_identity = percent_identity
        self.alignment_length = alignment_length
        self.mismatches = mismatches
        self.gap_opens = gap_opens
        self.query_start = query_start
        self.query_end = query_end
        self.subject_start = subject_start
        self.subject_end = subject_end
        self.evalue = evalue
        self.bit_score = bit_score

    @classmethod
    def parse_blast_list(cls, l):
        return cls(
            *(type_(v) for (type_, v) in zip(cls._field_types, l))
        )

    @classmethod
    def parse_blast_line(cls, l):
        return cls.parse_blast_list(l.split())

    @classmethod
    def parse_blast_lines(cls, lines):
        return [cls.parse_blast_line(l) for l in lines if not l.startswith("#")]
                
                
