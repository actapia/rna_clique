import functools

import numpy as np
import pandas as pd

from functools import cached_property
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Optional, Any
from numbers import Real
from fractions import Fraction

from gene_matches_tables import read_table
from multiset_key_dict import MultisetKeyDict
from identity import id_

def similarities_from_dfs(
    tables: Iterable[tuple[frozenset[str], pd.DataFrame]]
) -> Iterable[tuple[frozenset[str], Fraction]]:
    """Yield similarities for pairs of samples using gene matches tables.

    Each value yielded is a pair. The first element of the pair is itself a
    frozenset containing the IDs of the two samples for which the similarity was
    computed. The second element is a Fraction object, the similarity between
    the two samples.

    This function can raise a ZeroDivsionError when attempting to yield the
    similarity for a pair of samples if the filtered gene matches table for that
    sample pair is an empty dataframe. In that case, the similarity could be
    considered undefined or unknown.
    """
    
    for (qsample, ssample), restricted in tables:
        dist = Fraction(
            int(restricted["nident"].sum()),
            int(restricted["length"].sum() - restricted["gaps"].sum())
        )
        yield frozenset((qsample, ssample)), dist

class ComparisonSimilarityComputer:
    """Base class for computing similarities from comparison statistics.

    The comparion_dfs provide an implicit mapping between unordered pairs of
    samples to the dataframes containing tables of statistics comparing the
    samples (based, for example, on BLAST results). The "keys" in the mapping
    are expected to be frozensets with two string elements, the names of the
    pair of samples. The corresponding value for a given pair of samples is the
    Pandas dataframe containing per-sequence comparison (e.g., alignment)
    statistics for that pair of samples. A MultisetKeyDict object may be
    provided as the comparison_dfs, or the generator returned by the items
    method of an ordinary dict may be used instead.

    This class is not be used directly---use one of its subclasses instead.

    Attributes:
        comparison_dfs: An iterable mapping sample pairs to comparisons.
    """
    
    # Columns of comparison_dfs for which to use categorical data types.
    categorical_columns = []
    
    def __init__(
            self,
            comparison_dfs: Iterable[tuple[frozenset[str, str], pd.DataFrame]],
            sample_count: Optional[int] = None
    ):
        """Construct a ComparisonSimilarityComputer for the given comparisons.

        The number of samples can be provided explicitly. If it is not, the
        sample_count property will lazily determine the number of samples from
        the comparisons when needed.

        This class does not provide implementations for all necessary methods,
        so it should not be instantiated directly. Instead, use one of the
        classes that subclass this one.

        Parameters:
            comparison_dfs:     Mapping from sample pairs to comparison tables.
            sample_count (int): Number of samples in the analysis.
        """
        self.comparison_dfs = comparison_dfs
        self._sample_count =  sample_count
        self._samples = None

    @classmethod
    def mapping_from_dfs(
            cls,
            dfs: Iterable[pd.DataFrame]
    ) -> Iterator[tuple[frozenset[str], pd.DataFrame]]:
        """Iterate over pairs of samples and their comparison tables.

        This function takes an Iterable of dataframes and yields for each
        dataframe a (ordered) pair consisting of the unordered pair of samples
        for which the dataframe provides comparison statistics and the dataframe
        itself. This provides an implicit mapping from unordered pairs of
        samples to dataframes; the mapping can be made explicit by constructing
        a dict or similar Mapping using the yielded key-value pairs.

        This function obtains each dataframe's samples by examining the
        "qsample" and "ssample" columns of the dataframes. These columns are
        assumed to be constant for each dataframe, so the sample names are taken
        to be the values in the columns' first rows.

        Parameters:
            dfs: Iterable of comparison dataframes for different sample pairs.
        """
        for df in dfs:
            qsample = df["qsample"][0]
            ssample = df["ssample"][0]
            yield frozenset((qsample, ssample)), df

    @classmethod
    def _read_table(
            cls,
            table_path: Path,
            remove_seqids: bool = True,
            convert_to_categorical: bool = True,
    ) -> pd.DataFrame:
        """Read a comparison table from a file using appropriate options.

        This function behaves similarly to the gene_matches_tables.read_table
        function but includes some additional processing appropriate for
        comparison tables used by instances of this class.

        This function can take some steps to reduce the memory usage of the
        loaded dataframe. First, it can remove the qseqid and sseqid columns,
        since they are typically redundant with the parsed qgene, qiso, sgene,
        and siso columns, which, being integers, use less space. This behavior
        can be disabled by passing False for the remove_seqids parameter.

        The function can also convert certain columns to categorical
        datatypes. The columns to be converted are specified in the
        categorical_columns attribute of the class itself. This behavior can be
        disabled by passing False for the convert_to_categoricl argument.

        Parameters:
            table_path:                    Path to comparisons table to load.
            remove_seqids (bool):          Remove qseqid and sseqid columns.
            convert_to_categorical (bool): Convert certain columns to
                                           categorical datatype.

        Returns:
            The table at the specified path, with needed changes applied.
        """
        table = read_table(table_path)
        if remove_seqids:
            for col in ["q", "s"]:
                try:
                    table = table.drop(col + "seqid", axis=1)
                except KeyError:
                    pass
        if convert_to_categorical:
            cols = [c for c in cls.categorical_columns if c in table.index]
            table[
                cols
            ] = table[
                cols
            ].astype("category")
        return table
            
    @classmethod
    def _load_tables(
            cls,
            table_paths: Iterable[Path],
            *args,
            **kwargs
    ) -> Iterator[tuple[frozenset[str, str], pd.DataFrame]]:
        """Load multiple tables with given settings and return implicit mapping.

        This function uses the _read_table function of the same class to read
        multiple comparison tables and produces an iterator over pairs of
        samples and their corresponding loaded comparison dataframes using the
        class's mapping_from_dfs function.

        Parameters:
            table_paths: Paths to comparison tables to load.
        """
        return cls.mapping_from_dfs(
            map(
                functools.partial(cls._read_table, *args, **kwargs),
                table_paths
            )
        )

    @property
    def sample_count(self):
        """The number of samples in the similarity matrix."""
        if self._sample_count is None:
            self._sample_count = len(
                set.union(
                    *(
                        {a, b} for ((a, b), df) in self.comparison_dfs
                    )
                )
            )
        return self._sample_count

    def _similarity_helper(self) -> Iterator[tuple[frozenset[str], Real]]:
        """Get implicit mapping from sample pairs to similarities.

        This method should return a generator yielding (ordered) pairs. The
        first element of each pair should be an unordered pair (i.e., a
        frozenset) of strings representing the samples for which the similarity
        should be computed. The second element should be the similarity itself.

        This class has no implementation for _similarity_helper---see the
        subclasses instead.
        """
        raise NotImplementedError()

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

    @classmethod
    def similarity_to_dissimilarity(cls, sim: Real) -> Real:
        """Obtain a dissimilarity (or distasnce) from a similarity.

        By default, this value is just 1 minus the similarity.

        Parameters:
            sim: The similarity for which to get corresponding dissimilarity.

        Returns:
            The dissimilarity (distance) corresponding to the similarity value.
        """
        return 1 - sim

    def get_dissimilarities(self) -> MultisetKeyDict[Any, Real]:
        """Returns the dissimilarities between pairs of samples."""
        #res = MultisetKeyDict(self.similarities)
        #for k, v in res.
        return MultisetKeyDict(
            {
                p: self.similarity_to_dissimilarity(s)
                for (p, s) in self.similarities.items()
            }
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
        """Make a matrix from a mapping from unordered pairs to values.

        Each row and column in the result matrix corresponds to a sample. The
        order of the samples in the rows and columns is exactly their order in
        the samples property. The element at the row corresponding to sample A
        and column corresponding to sample B is the entry in the provided dict
        for that pair of samples.

        Parameters:
            d: A MultisetKeyDict mapping unordered sample pairs to real values.

        Returns:
            An equivalent matrix for the provided mapping.
        """        
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
        """Convert a matrix of values for pairs of samples to a dataframe.

        This method assumes that the rows and columns in the provided matrix
        both correspond to the samples in the samples property, in the same
        order.

        Parameters:
            mat: The numpy array to conver to a dataframe.

        Returns:
            An equivalent Pandas dataframe, with labeled rows and columns.
        """
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
    def _constructor_args_from_filenames(
            cls,
            comparison_fns: Iterable[Path],
            store_dfs: bool = True,
            *args,
            **kwargs
    ) -> tuple[list, dict[str, Any]]:
        """Get constructor arguments for constructing from filenames.

        Parameters:
            comparison_fns:   Iterable of Paths to gene matches tables.
            store_dfs (bool): Whether to store gene matches tables in memory.

        Returns:
            The positional and keyword constructor arguments.
        """
        if store_dfs:
            f = MultisetKeyDict
        else:
            f = id_
        args = [f(cls._load_tables(comparison_fns))] + list(args)
        return args, kwargs

    @classmethod
    def from_filenames(
            cls,
            *args,
            **kwargs
    ):
        """Constructs a SimilarityComputer from paths to table files.

        To save memory, the store_dfs parameter can be set to False. In that
        case, comparison_dfs will be generator, and it will not be possible to
        access them more than once.

        Parameters:
            comparison_fns:   Paths to stored gene matches tables.
            store_dfs (bool): Whether to store the dataframes loaded.

        Returns:
            A SimilarityComputer using the given gene matches graph and tables.
        """
        args, kwargs = cls._constructor_args_from_filenames(*args, **kwargs)
        return cls(*args, **kwargs)

