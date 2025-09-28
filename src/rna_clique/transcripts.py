import re

from typing import Callable, Any
from collections import namedtuple

# This default regex is based on the transcript ID format used in rnaSPAdes
# 4.0.0.
default_gene_re = re.compile(r"^.*cov_([0-9]+(?:\.[0-9]+))_g([0-9]+)_i([0-9]+)")

def cast_namedtuple(
        name: str,
        field_names: list[str],
        converters: list[Callable[[Any], Any]],
        **kwargs
) -> type[tuple]:
    """Get a namedtuple-like class that converts its attributes on construction.

    This function behaves similarly to namedtuple and ultimately uses that
    function, but this function also gives the created class a custom __new__
    function that uses the provided converters parameter to "cast" the arguments
    to the correct types when constructing a new instance of the class.

    Like namedtuple, this function requires the name of the class to create and
    the names of the attributes, in order. This function additionally requires
    the converts parameter, a parallel list containing the unary functions to
    apply to constructor arguments.

    Additional arbitrary keyword arguments that can be passed to namedtuple are
    also supported.

    Parameters:
        name (str):         Name of the namedtuple class to create.
        field_names (list): Names of fields of the tuple, in order.
        converters (list):  Functions to convert field values in constructor.

    Returns:
        A namedtuple that converts attributes using provided converters.
    """
    class _base(namedtuple(name + "_base", field_names, **kwargs)):
        def __new__(cls, *args, **kwargs):
            kwargs |= dict(zip(cls._fields, args))
            return super().__new__(
                cls,
                **{f: c(kwargs[f]) for (f, c) in zip(cls._fields, converters)}
            )
    return type(name, (_base,), {})

# TranscriptID is a namedtuple containing metadata about a transcript.
#
# The coverage attribute is the k-mer coverage of the transcript, a float value.
# It quantifies how much of the input RNA-seq data contributes to the
# transcript.
#
# The gene attribute is an integer identifying the gene to which the transcript
# belongs.
#
# The isoform attribute is an integer identifying which isoform of the gene the
# transcript is.
TranscriptID = cast_namedtuple(
    "TranscriptID",
    ["coverage", "gene", "isoform"],
    [float, int, int]
)

def re_parse_transcript_id(
        cls: type,
        expr: re.Pattern
) -> Callable[[str], TranscriptID]:
    """Create a function that parses transcript FASTA IDs using a regex.

    To offer flexibility, the provided regular expression can use either
    positional or named capture groups for the TranscriptID fields. When
    positional capture groups are used, they are assumed to correspond to
    coverage, gene, and isoform, in that order.

    A combination of positional and named capture groups can also be used. In
    that case, the positional groups are assumed to be the non-named groups, in
    the same order above.

    Parameters:
        expr: Regular expression for parsing transcript FASTA IDs.

    Returns:
        A function that uses the regex to parse a FASTA ID into a TranscriptID.
    """
    def parse_transcript_id(id_: str) -> TranscriptID:
        """Parse the given FASTA ID into a TranscriptID object using a regex.

        This function was created using the re_parse_transcript_id function; the
        regex this function uses was provided to that function.

        Parameters:
            id_ (str): The FASTA ID to parse.

        Returns:
            The TranscriptID object parsed from the given FASTA ID.
        """
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

default_parser = TranscriptID.parser_from_re(default_gene_re)
