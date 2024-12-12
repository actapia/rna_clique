import argparse
import pickle
import sys
import functools

import numpy as np
import pandas as pd
import networkx as nx

from functools import cached_property
from fractions import Fraction
from pathlib import Path
from typing import Optional, Any
from collections.abc import Iterable, Iterator
from numbers import Real

from build_graph import component_subgraphs
from multiset_key_dict import MultisetKeyDict, FrozenMultiset
from gene_matches_tables import read_table

from tqdm import tqdm

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

def id_(x):
    """Identity function."""
    return x

def handle_arguments():
    parser = argparse.ArgumentParser(
        description="compute pairwise similarities from graph and comparisons"
    )
    parser.add_argument(
        "-g",
        "--graph",
        type=Path,
        required=True,
        help="path to the gene matches graph pickle"
    )
    parser.add_argument(
        "-c",
        "--comparisons",
        type=Path,
        nargs="+",
        required=True,
        help="paths to the gene matches tables"
    )
    parser.add_argument(
        "-s",
        "--samples",
        type=int,
        help="number of samples"
    )
    parser.add_argument(
        "-e",
        "--embed",
        action="store_true",
        help="enter an IPython shell after computing the matrix"
    )
    parser.add_argument(
        "-o",
        "--out-type",
        choices=["sim", "dis"], # Similarity, dissimilarity
        default="sim",
        help="type of matrix to produce (similarity or dissimilarity)"
    )
    parser.add_argument(
        "-l",
        "--print-sample-list",
        action="store_true",
        help="print the list of samples before the matrix"
    )
    return parser.parse_args()

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

def print_mat(m : np.ndarray):
    np.savetxt(sys.stdout, m, fmt="%s", delimiter=' ')



    
class SampleSimilarity:
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

    The comparion_dfs provide an implicit mapping between unordered pairs of
    samples to the dataframes containing the BLAST results of the comparisons
    between the samples. The "keys" in the mapping are expected to be frozensets
    with two string elements, the names of the pair of samples. The
    corresponding value for a given pair of samples is the Pandas dataframe
    containing BLAST results for that pair of samples. A MultisetKeyDict object
    may be provided as the comparison_dfs, or the generator returned by the
    items method of an ordinary dict may be used instead.

    If the graph is saved in a pickle, and the comparison dataframes are stored
    in pickles or HDF files, the from_filenames classmethod may provide a more
    convenient way of constructing a SampleSimilarity object.

    Attributes:
        graph:         The gene matches graph representing gene orthologies.
        comparison_dfs: An iterable mapping sample pairs to comparisons.
    """
    def __init__(
            self,
            graph: nx.Graph,
            comparison_dfs: Iterable[tuple[frozenset[str, str], pd.DataFrame]],
            sample_count: Optional[int] = None
    ):
        self.graph = graph
        self.comparison_dfs = comparison_dfs
        self._sample_count = sample_count
        self._samples = None

    @classmethod
    def mapping_from_dfs(cls, dfs):
        for df in dfs:
            qsample = df["qsample"][0]
            ssample = df["ssample"][0]
            yield frozenset((qsample, ssample)), df
            

    @classmethod
    def _load_tables(
            cls,
            table_paths: Iterable[Path]
    ) -> Iterator[tuple[frozenset[str, str], pd.DataFrame]]:
        return cls.mapping_from_dfs(map(read_table, table_paths))
            

    @property
    def sample_count(self):
        """The number of samples in the similarity matrix."""
        if self._sample_count is None:
            self._sample_count = len(
                set(
                    n[0] for comp in component_subgraphs(
                        self.graph
                    ) for n in comp.nodes
                )
            )
        return self._sample_count

    @cached_property
    def valid(self):
        """A dataframe containing all genes found in ideal components."""
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
        return functools.reduce(
            functools.partial(
                restrict_to,
                self.valid
            ),
            [
                [a + b for b in ["sample", "gene"]]
                for a in ["s", "q"]
            ],
            comp_df
        )

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

    def _similarity_helper(self):
        for (qsample, ssample), comp_df in self.comparison_dfs:
            restricted = self.restricted(comp_df)
            dist = Fraction(
                restricted["nident"].sum(),
                restricted["length"].sum() - restricted["gaps"].sum()
            )
            yield frozenset((qsample, ssample)), dist

    @cached_property
    def similarities(self):
        """The similarities between pairs of samples."""
        sims = MultisetKeyDict(self._similarity_helper())
        self._samples = sorted(list(sims.key_elements()))
        #print(sims._dict)
        res = sims | MultisetKeyDict(
            {(s, s): 1 for s in self._samples}
        )
        #print(res._dict)
        return res

    def get_dissimilarities(self) -> MultisetKeyDict[Any, Real]:
        """Returns the dissimilarities between pairs of samples."""
        #res = MultisetKeyDict(self.similarities)
        #for k, v in res.
        return MultisetKeyDict(
            {p: 1 - s for (p, s) in self.similarities.items()}
        )

    def get_similarities(self) -> MultisetKeyDict[Any, Real]:
        """Returns the similarities between pairs of samples."""
        return self.similarities

    @property
    def samples(self):
        """The samples in the analysis.

        Since the list of samples is obtained in the process of computing
        similarities, this property will compute similarities first if they have
        not been computed already. If the similarities are not needed, this
        property is an inefficient way of obtaining the list of samples.
        """
        if self._samples is None:
            self.similarities
        return self._samples

    def _pair_dict_to_matrix(
            self,
            d : MultisetKeyDict[Any, Real]
    ) -> np.ndarray:
        return np.vstack(
            [
                [
                    float(
                        d[[a, b]]
                    ) for b in self.samples
                ] for a in self.samples
            ]
        )

    def get_similarity_matrix(self) -> np.ndarray:
        """Returns the computed pairwise similarity matrix.

        The order of the rows and columns is given by the samples property.

        Returns:
            A matrix giving the similarity for each pair of samples.
        """
        return self._pair_dict_to_matrix(self.get_similarities())

    def _matrix_to_df(self, mat: np.ndarray) -> pd.DataFrame:
        return pd.DataFrame(
            mat
        ).set_axis(
            self.samples,
            axis=1
        ).set_axis(
            self.samples,
            axis=0
        )

    def get_similarity_df(self) -> pd.DataFrame:
        """Returns the pairwise similarity matrix as a Pandas dataframe.

        Rows and columns are indexed with the sample names.

        Returns:
            A Pandas dataframe giving pairwise similarities for all samples.
        """
        return self._matrix_to_df(self.get_similarity_matrix())

    def get_dissimilarity_matrix(self) -> np.ndarray:
        """Returns the computed pairwise dissimilarity (distance) matrix.

        The order of the rows and columns is given by the samples property.

        Returns:
            A matrix giving the dissimilarity for each pair of samples.
        """
        return self._pair_dict_to_matrix(self.get_dissimilarities())

    def get_dissimilarity_df(self) -> pd.DataFrame:
        """Returns the pairwise dissimilarity matrix as a Pandas dataframe.

        Rows and columns are indexed with the sample names.

        Returns:
            A Pandas dataframe giving pairwise dissimilarities for all samples.
        """
        return self._matrix_to_df(self.get_dissimilarity_matrix())

    @classmethod
    def from_filenames(
            cls,
            graph_fn : Path,
            comparison_fns : Iterable[Path],
            store_dfs : bool = True,
            *args,
            **kwargs
    ):
        """Constructs a SampleSimilarity from paths to graph and table files.

        Parameters:
            graph_fn:         Path to pickle for gene matches graph.
            comparison_fns:   Paths to stored gene matches tables.
            store_dfs (bool): Whether to store the dataframes loaded.

        Returns:
            A SampleSimilarity using the pickled gene matches graph and tables.
        """
        with open(graph_fn, "rb") as f:
            graph = pickle.load(f)
        if store_dfs:
            f = MultisetKeyDict
        else:
            f = id_
        #embed()
        return cls(graph, f(cls._load_tables(comparison_fns)), *args, **kwargs)
    
        
def main():
    args = handle_arguments()
    sim = SampleSimilarity.from_filenames(
        args.graph,
        tqdm(args.comparisons),
        sample_count=args.samples
    )
    if args.print_sample_list:
        print("\n".join(sim.samples))
    if args.out_type == "sim":
        print_mat(sim.get_similarity_matrix())
    elif args.out_type == "dis":
        print_mat(sim.get_dissimilarity_matrix())
    if args.embed:
        from IPython import embed
        embed()

        
if __name__ == "__main__":
    main()
