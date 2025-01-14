from fractions import Fraction
from similarity_computer import ComparisonSimilarityComputer

class UnfilteredSimilarity(ComparisonSimilarityComputer):
    def _similarity_helper(self):
        for (qsample, ssample), comp_df in self.comparison_dfs:
            dist = Fraction(
                comp_df["nident"].sum(),
                comp_df["length"].sum() - comp_df["gaps"].sum()
            )
            yield frozenset((qsample, ssample)), dist
