import re
from collections import namedtuple

default_gene_re = re.compile(r"^.*cov_([0-9]+(?:\.[0-9]+))_g([0-9]+)_i([0-9]+)")

def cast_namedtuple(name, field_names, converters, **kwargs):
    class _base(namedtuple(name + "_base", field_names, **kwargs)):
        def __new__(cls, *args, **kwargs):
            kwargs |= dict(zip(cls._fields, args))
            return super().__new__(
                cls,
                **{f: c(kwargs[f]) for (f, c) in zip(cls._fields, converters)}
            )
    return type(name, (_base,), {})

TranscriptID = cast_namedtuple(
    "TranscriptID",
    ["coverage", "gene", "isoform"],
    [float, int, int]
)

def re_parse_transcript_id(cls: type, expr: re.Pattern):
    def parse_transcript_id(id_: str):
        d = {}
        remaining = []
        match_ = expr.search(id_)
        for field in cls._fields:
            try:
                d[field] = match_.group(field)
            except IndexError:
                remaining.append(field)
        named = set(expr.groupindex.values())
        for field, i in zip(
                remaining,
                (i for i in range(1, expr.groups + 2) if i not in named)
        ):
            try:
                d[field] = match_.group(i)
            except IndexError as e:
                print("Bad ID?", id_)
                raise e
        return TranscriptID(**d)
    return parse_transcript_id

TranscriptID.parser_from_re = classmethod(re_parse_transcript_id)
