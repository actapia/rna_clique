import argparse
import sys
import pickle
import re
import functools
from simple_blast.blasting import BlastnSearch
import numpy as np

from fractions import Fraction

import pandas as pd

from typing import Callable

from IPython import embed

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("transcripts1")
    parser.add_argument("transcripts2")
    parser.add_argument(
        "--regex",
        "-r",
        type=re.compile,
        default=re.compile("^.*g([0-9]+)_i([0-9]+)")
    )
    parser.add_argument("-e", "--evalue", type=float, default=1e-50)
    parser.add_argument(
        "-n",
        "--top-n",
        type=int,
        default=1,
        help="use top n matches"
    )
    parser.add_argument(
        "--keep-all",
        "-k",
        action="store_true",
        help="keep all pairs in case of a tie"
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="hide presumed ortholog distances"
    )
    parser.add_argument(
        "-f",
        "--report-float",
        action="store_true",
        help="report float instead of fraction"
    )
    return parser.parse_args()

def parse_seq_id(regex: re.compile, s: pd.Series) -> pd.DataFrame:
    """Parse a seq_id column to extract the gene and isoform IDs."""
    return s.str.extract(regex).astype(np.int32)

def gene_matches(
        parse: Callable[[re.compile, pd.Series], pd.DataFrame],
        path1: str,
        path2: str,
        evalue: float,
        n: int = 1,
        **blast_kwargs
) -> pd.DataFrame:
    """Find the isotigs in file 2 that best match each gene in file 1.

    This function performs a BLAST search of the sequences in the FASTA file
    located at path1 again the sequences in the FASTA file located at path2.

    For each gene (isotig set) g_1_x in the first file, the top n pairs of 
    isotigs (i_1_a, i_2_b) such that i_1_a in g_1_x is found, where pairs are
    compared by their bitscores in the BLAST search. The genes to which the
    second elements of these pairs belongs may be treated as candidate homologs
    of g_1_x in sample 2.

    This function must be provided a function to parse the sequence IDs and
    extract the numerical gene and isotig IDs.

    This function accepts an evalue to be used as a cutoff in the BLAST search.
    Additional parameters for the BLAST search may also be provided as variadic
    arguments to this function.

    Parameters:
        parse:          Function to parse sequence IDs into gene and isotig IDs.
        path1 (str):    Path to the FASTA file used as the BLAST search query.
        path2 (str):    Path to the FASTA file used as the BLAST search subject.
        evalue (float): Expect value cutoff for BLAST search.
        n (int):        Number of top isotigs to select for each query gene.

    Returns:
        A dataframe of BLAST hits for subject isotigs best matching query genes.
    """
    # TODO: Check whether we really want the top n subject isotigs or the top
    # n subject genes. (I suspect we really want the latter.)
    search = BlastnSearch(path1, path2, evalue=evalue, **blast_kwargs)
    for t in ["q", "s"]:
        search.hits[[t + "gene", t + "iso"]] = parse(search.hits[t + "seqid"])
    return highest_bitscores(search.hits, n, keep="all")

eprint = functools.partial(print, file=sys.stderr)

def highest_bitscores(
        df: pd.DataFrame,
        n: int = 1,
        groupby: str = "qgene",
        **kwargs
) -> pd.DataFrame:
    """Select the rows (hits) with the highest bitscore for each group.
    
    This function accepts a dataframe containing BLAST hits with bitscores.
    By default, this function selects the top hit for each query gene, but
    the column on which to group and the number of hits to select per group may
    be controlled with the optional parameters.

    Additional keyword arguments will be passed to the pandas nlargest function
    to control its behavior. This allows the caller to change the strategy for
    breaking ties, for example.

    Parameters:
        df:            A dataframe containing BLAST hits.
        n (int):       Number of top hits to select per group.
        groupby (str): Column on which to group hits.
    """
    return df.loc[
        df.groupby(
            groupby
        )["bitscore"].nlargest(n, **kwargs).index.get_level_values(-1)
    ]

class HomologFinder:
    merge_columns = ["qgene", "sgene"]
    
    def __init__(
            self,
            regex: str,
            top_n: int,
            evalue: float,
            keep_all: bool,
            debug: bool = False,
            **blast_kwargs
    ):
        # self.regex = regex
        # self.top_n = top_n
        # self.evalue = evalue
        # Partially apply some functions to make the code below less repetitive.
        parse = functools.partial(parse_seq_id, regex)
        self.gm = functools.partial(
            gene_matches,
            parse=parse,
            evalue=evalue,
            n=top_n,
            additional_columns=["gaps", "nident"],
            **blast_kwargs
        )
        self.keep_all = keep_all
        self.debug = debug

    def get_match_table(
            self,
            transcripts1: str,
            transcripts2: str
    ) -> pd.DataFrame:
        # Search in both the forward and reverse direction to get the isotigs
        # (and corresponding gene IDs) in one sample that best match genes in
        # the other.
        if self.debug:
            eprint("Getting forward matches.")
        forward_matches = self.gm(
            path1=transcripts1,
            path2=transcripts2
        )
        if self.debug:
            eprint("Getting backward matches.")
        backward_matches = self.gm(
            path1=transcripts2,
            path2=transcripts1
        )
        # We rename the columns in the reverse matches to enable merging.
        backward_matches.rename(
            columns={
                a+v : b+v
                for (a,b) in [("q","s"),("s","q")]
                for v in ["seqid","gene","iso"]
            },
            inplace=True
        )
        # Compute the "intersection" of the forward and reverse matches.
        #
        # A row (hit) in either dataframe with query gene ID q and subject gene ID s
        # is kept if and only if there is a row with query gene ID q and subject
        # gene ID s in both dataframes.
        #
        # The resulting dataframe has two columns---index_x and index_y---that
        # indicate where in the original dataframes the rows were found.
        #
        # (Keep in mind that we swapped the order order of query and subject in the
        # reverse dataframe, so "query" is always something in sample 1, and
        # "subject" is always something in sample 2 from now on.)
        if forward_matches.empty or backward_matches.empty:
            intersection = pd.DataFrame(columns=self.merge_columns + ["index_x", "index_y"])
        else:
            intersection = pd.merge(
                forward_matches[self.merge_columns].reset_index(),
                backward_matches[self.merge_columns].reset_index(),
                how="inner",
                on=self.merge_columns
            )
        # Get the hits (and subject gene ID) with the highest bitscore for each
        # query gene.
        return highest_bitscores(
            # Get the hits with highest bitscore for each pair of homologous genes.
            highest_bitscores(
                # Get the original rows back from the forward and reverse
                # dataframes.
                pd.concat(
                    [
                        forward_matches.loc[
                            intersection["index_x"].drop_duplicates()
                        ],
                        backward_matches.loc[
                            intersection["index_y"].drop_duplicates()
                        ],
                    ]
                ).reset_index(drop=True),
                groupby=self.merge_columns,
                # TODO: Decide if "all" is the right choice here.
                keep="all"
            ),
            keep=["first", "all"][self.keep_all]
        )        
        # Compute distances at a gene level.
        # best_matches["dist"] = \
        #     best_matches["nident"] / \
        #     (best_matches["length"] - best_matches["gaps"])
        # Ignore rows with identical query and subject gene IDs.
        # return best_matches[merge_columns + ["dist"]].drop_duplicates()

    @classmethod
    def without_duplicates(
            cls,
            df: pd.DataFrame,
            columns: list[str] = None
    ) -> pd.DataFrame:
        if columns is None:
            columns = cls.merge_columns
        return df[columns].drop_duplicates()

def main():
    args = parse_arguments()
    match_finder = HomologFinder(
        args.regex,
        args.top_n,
        args.evalue,
        args.keep_all
    )
    best_matches = match_finder.get_match_table(
        args.transcripts1,
        args.transcripts2,        
    )
    dedup = match_finder.without_duplicates(best_matches)
    if not args.quiet:
        for match in dedup.itertuples(index=False):
            print(*match)
    eprint(f"Found {len(dedup)} matches.")
    # Compute the overall distance using all identified homologs.
    #
    # I use a Fraction here to avoid loss of precision.
    dist = Fraction(
        best_matches["nident"].sum(),
        best_matches["length"].sum() - best_matches["gaps"].sum()
    )
    if args.report_float:
        print(float(dist))
    else:
        print(dist)
    #embed()
    
if __name__ == "__main__":
    main()
