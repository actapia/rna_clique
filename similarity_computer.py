import numpy as np
import pandas as pd

from functools import cached_property
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Optional, Any
from numbers import Real

from gene_matches_tables import read_table
from multiset_key_dict import MultisetKeyDict

class ComparisonSimilarityComputer:
    """Base class for computing similarities from BLAST comparsions.

    The comparion_dfs provide an implicit mapping between unordered pairs of
    samples to the dataframes containing the BLAST results of the comparisons
    between the samples. The "keys" in the mapping are expected to be frozensets
    with two string elements, the names of the pair of samples. The
    corresponding value for a given pair of samples is the Pandas dataframe
    containing BLAST results for that pair of samples. A MultisetKeyDict object
    may be provided as the comparison_dfs, or the generator returned by the
    items method of an ordinary dict may be used instead.

    Attributes:
        comparison_dfs: An iterable mapping sample pairs to comparisons.
    """
    def __init__(
            self,
            comparison_dfs: Iterable[tuple[frozenset[str, str], pd.DataFrame]],
            sample_count: Optional[int] = None
    ):
        self.comparison_dfs = comparison_dfs
        self._sample_count =  sample_count
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
                set.union(
                    *(
                        {a, b} for ((a, b), df) in self.comparison_dfs
                    )
                )
            )
        return self._sample_count

    def _similarity_helper(self):
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
