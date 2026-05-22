import pickle
import functools
import sys

import pandas as pd
import networkx as nx

import config as config_module

from functools import cached_property
from fractions import Fraction
from pathlib import Path
from typing import Optional, Any
from collections.abc import Iterable, Iterator

from graph import component_subgraphs
from multiset_key_dict import FrozenMultiset
from gene_matches_tables import get_table_files
from similarity_computer import (
    ComparisonSimilarityComputer,
    similarities_from_dfs
)
from app import eprint, set_except_hook

def is_complete(g : nx.Graph) -> bool:
    """Returns whether g is a complete (sub)graph."""
    v = len(g)
    return 2*len(g.edges) == v*(v-1)

def get_ideal_components(
        g : nx.Graph,
        samples: int
) -> Iterator[nx.Graph]:
    """Yields the ideal components of g, assuming a given number of samples."""
    for s in component_subgraphs(g):
        if len(s) == samples and is_complete(s):
            yield s

def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager(
        description=(
            "Compute pairwise distances from gene matches tables and graph."
        ),
    )
    arg_config.expose_fields_with_default_aliases(
        "graph",
        "tables_dir",
        "matrix",        
        required=True
    )
    arg_config.expose_fields_with_default_aliases(
        "output_dir",
    )

    # arg_config.add_argument(
    #     "-e",
    #     "--embed",
    #     action="store_true",
    #     help="enter an IPython shell after computing the matrix"
    # )
    arg_config.add_output_config_argument()
    return arg_config

def restrict_to(
        df2 : pd.DataFrame,
        df1 : pd.DataFrame,
        columns : Iterable[str]
) -> pd.DataFrame:
    """Returns the rows of df1 for which the column values appear in df2.

    This function returns a "restricted" version of df1 in which only rows that
    have values for the provided columns which appear in df2 are kept.

    In other words, for each row r of df1, the values for the specified columns
    (r[columns]) are considered. If there is some row s in d2 such that
    tuple(r[columns]) == tuple(s), then the row is kept. 

    The columns provided should be a subset of the columns in df1 and correspond
    exactly to the columns in df2, though the columns in df1 and df2 need not
    have the same names.

    Parameters:
        df2:     The dataframe used for filtering.
        df1:     The dataframe from which to draw rows.
        columns: The columns of df1 that correspond to those of df2.

    Returns:
        df1 restricted to rows for which the values for columns appear in df2.
    """
    return df1.loc[
        pd.merge(
            df1[columns].rename(
                columns=dict(zip(columns, df2.columns))
            ).reset_index(),
            df2.reset_index(),
            how="inner",
            on=list(df2.columns)
        )["index_x"]
    ]

# TODO: Clarify the explanation in this docstring.
def restrict_multi(
        df2: pd.DataFrame,
        df1: pd.DataFrame,
        columns: Iterable[Iterable[str]]
):
    """Restrict df1 to rows in df2 based on multiple lists of columns.

    This function essentially performs restrict_to multiple times with the same
    df2 but different columns in each step. The dataframe to be filtered is
    initially df1, and the restricted dataframe after each step is used as the
    dataframe to be filtered for the next.

    Parameters:
        df2:     The dataframe used for filtering.
        df1:     The dataframe from which to draw rows.
        columns: The lists of columns of df1 that correspond to those of df2.

    Returns:
        df1, without rows where values for some column list isn't a row in df2.
    """
    return functools.reduce(functools.partial(restrict_to, df2), columns, df1)

class NoIdealComponentsError(Exception):
    pass
    
class SampleSimilarity(ComparisonSimilarityComputer):
    """Computes samples' similarities from gene matches graph and comparisons.

    In the gene matches graph, each vertex represents a gene in a specific
    sample. Hence, SampleSimilarity expects each vertex to be a pair. The first
    element of the pair should be a string representing the sample. The name
    used should be the same as those used in the qsample and ssample columns of
    the comparison dataframes. The second element of the pair should be an
    integer representing the gene ID and should likewise correspond to the
    values for the qgene or sgene columns in the comparison dataframes.

    Edges in the gene matches graph should exist between pairs of genes inferred
    to be orthologs.

    If the graph is saved in a pickle, and the comparison dataframes are stored
    in pickles or HDF files, the from_filenames classmethod may provide a more
    convenient way of constructing a SampleSimilarity object.

    Attributes:
        graph:          The gene matches graph representing gene orthologies.
        comparison_dfs: An iterable mapping sample pairs to comparisons.
    """

    # List of lists of columns corresponding to sample and gene IDs for subject
    # and queries.
    sample_gene_columns = [
        [a + b for b in ["sample", "gene"]]
        for a in ["s", "q"]
    ]

    # Columns that can be stored as Pandas categorical values.
    categorical_columns = ["qsample", "ssample", "sstrand"]
    
    def __init__(
            self,
            graph: nx.Graph,
            comparison_dfs: Iterable[tuple[frozenset[str, str], pd.DataFrame]],
            sample_count: Optional[int] = None
    ):
        super().__init__(comparison_dfs, sample_count)
        self.graph = graph
            
    @property
    def sample_count(self):
        """The number of samples in the similarity matrix."""
        if self._sample_count is None:
            self._sample_count = len(self.samples)
        return self._sample_count

    @cached_property
    def valid(self):
        """A dataframe containing all genes found in ideal components."""
        #from IPython import embed; embed()
        return pd.DataFrame(
            (
                n for comp in get_ideal_components(
                    self.graph,
                    self.sample_count
                ) for n in comp.nodes
            ),
            columns=["sample", "gene"]
        )

    def restricted(self, comp_df: pd.DataFrame) -> pd.DataFrame:
        """Returns the provided dataframe, restricted to valid genes.

        Valid genes are those that are found in some ideal component.

        Parameters:
            comp_df: The dataframe to restrict to valid genes.

        Returns:
            comp_df, restricted to genes appearing in ideal components.
        """
        return restrict_multi(self.valid, comp_df, self.sample_gene_columns)

    def restricted_comparison_dfs(
            self
    ) -> Iterator[tuple[FrozenMultiset, pd.DataFrame]]:
        """Returns comparison_dfs, restricted to valid genes.

        Valid genes are those that are found in some ideal component.

        This function is a generator that yields key-value pairs that may be
        interpreted as an implicit mapping. To obtain a MultisetKeyDict for
        this mapping, simply construct a MultisetKeyDict from this generator.

        Each key is a FrozeznMultiset representing an unordered pair of
        samples. The value corresponding to a given pair of samples is the
        Pandas dataframe containing the BLAST results for that pair of samples,
        restricted to genes found in some ideal component.

        Returns:
            A generator mapping sample pairs to their restricted BLAST results.
        """       
        for k, df in self.comparison_dfs.multiset_iter():
            yield k, self.restricted(df)

    def _similarity_helper(self) -> Iterator[tuple[frozenset[str], Fraction]]:
        """Yield similarities for pairs of samples using filtered tables.

        This function can raise a NoIdealComponentsError when attempting to
        yield the similarity for a pair of samples if the filtered gene matches
        table for that sample pair is an empty dataframe. In that case, the
        similarity could be considered undefined or unknown.
        """
        try:
            yield from similarities_from_dfs(
                (k, self.restricted(v)) for (k, v) in self._comparison_df_iter
            )
        except ZeroDivisionError:
            raise NoIdealComponentsError()

    @classmethod
    def _constructor_args_from_filenames(
            cls,
            graph_fn : Path,
            comparison_fns : Iterable[Path],
            store_dfs : bool = True,
            remove_seqids: bool = True,
            convert_to_categorical: bool = True,
            *args,
            **kwargs
    ) -> tuple[list, dict[str, Any]]:
        """Get constructor arguments for constructing from filenames.

        To save memory, the store_dfs parameter can be set to False. In that
        case, comparison_dfs will be a generator, and it will not be possible to
        access them more than once.

        The qseqid and sseqid columns are often long and can also be removed to
        save memory. Likewise, certain columns (specified in the class's
        cateogrical_columns attribute) can be made Pandas categorical columns,
        which can further reduce the memory footprint.        

        Parameters:
            comparison_fns:                 Paths to stored gene matches tables.
            store_dfs (bool):               Store the dataframes loaded.
            remove_seqids (bool):           Delete seqid columns.
            convert_to_categorical (bool):  Make certain columns categorical.

        Returns:
            The positional and keyword constructor arguments.
        """        
        args, kwargs = super()._constructor_args_from_filenames(
            comparison_fns,
            store_dfs,
            *args,
            **kwargs
        )
        with open(graph_fn, "rb") as f:
            graph = pickle.load(f)
        args = [graph] + args
        return args, kwargs

    @classmethod
    def from_filenames(
            cls,
            *args,
            **kwargs
    ):
        """Constructs a SampleSimilarity from paths to table files.

        To save memory, the store_dfs parameter can be set to False. In that
        case, comparison_dfs will be generator, and it will not be possible to
        access them more than once.

        The qseqid and sseqid columns are often long and can also be removed to
        save memory. Likewise, certain columns (specified in the class's
        cateogrical_columns attribute) can be made Pandas categorical columns,
        which can further reduce the memory footprint.        

        Parameters:
            comparison_fns:                 Paths to stored gene matches tables.
            store_dfs (bool):               Store the dataframes loaded.
            remove_seqids (bool):           Delete seqid columns.
            convert_to_categorical (bool):  Make certain columns categorical.

        Returns:
            A SampleSimilarity using the given gene matches graph and tables.        
        """        
        return super().from_filenames(*args, **kwargs)
        
def main():
    with set_except_hook():
        _, args, config = build_parser().get_arguments_and_config()
    with set_except_hook(config.verbose):
        tables = list(get_table_files(config.tables_dir))
        if not tables:
            eprint(
                "Warning: No gene matches tables found in {}".format(
                    config.tables_dir
                )
            )        
        sim = SampleSimilarity.from_filenames(
            config.graph,
            tables,
        )
        try:
            mat = sim.get_dissimilarity_df()
            mat.to_hdf(config.matrix, key="matrix", mode="w")
            config.mark_finish()
        except NoIdealComponentsError:
            eprint("No ideal components found. Cannot report distances!")
            sys.exit(1)
        if args.output_config:
            config.yaml_save(args.output_config)

if __name__ == "__main__":
    main()
