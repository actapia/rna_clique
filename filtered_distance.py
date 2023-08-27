import argparse
import pickle
import sys
import numpy as np
import pandas as pd
import functools
from functools import cached_property
from fractions import Fraction
from pathlib import Path
from plot_component_sizes import component_subgraphs
from find_homologs import eprint
from IPython import embed
from tqdm import tqdm

def is_complete(g):
    v = len(g)
    return 2*len(g.edges) == v*(v-1)

def get_ideal_components(g, samples):
    for s in component_subgraphs(g):
        if len(s) == samples and is_complete(s):
            yield s

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--graph", type=Path, required=True)
    parser.add_argument(
        "-c",
        "--comparisons",
        type=Path,
        nargs="+",
        required=True
    )
    parser.add_argument(
        "-s",
        "--samples",
        type=int
    )
    parser.add_argument(
        "-e",
        "--embed",
        action="store_true"
    )
    parser.add_argument(
        "-o",
        "--out-type",
        choices=["sim", "dis"], # Similarity, dissimilarity
        default="sim"
    )
    return parser.parse_args()

def restrict_to(df2, df1, columns):
    try:
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
    except KeyError as e:
        embed()

def print_mat(m):
    np.savetxt(sys.stdout, m, fmt="%s", delimiter=' ')

class SampleSimilarity:
    def __init__(self, graph, comparison_dfs, sample_count=None):
        self.graph = graph
        self.comparison_dfs = comparison_dfs
        self._sample_count = sample_count
        self._samples = None

    @classmethod
    def load_pickles(cls, pickles):
        for pick in pickles:
            yield pd.read_pickle(pick)

    @property
    def sample_count(self):
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
        return pd.DataFrame(
            (
                n for comp in get_ideal_components(
                    self.graph,
                    self.sample_count
                ) for n in comp.nodes
            ),
            columns=["sample", "gene"]
        )

    def _similarity_helper(self):
        for comp_df in self.comparison_dfs:
            qsample = comp_df["qsample"][0]
            ssample = comp_df["ssample"][0]
            restricted = functools.reduce(
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
            dist = Fraction(
                restricted["nident"].sum(),
                restricted["length"].sum() - restricted["gaps"].sum()
            )
            yield frozenset((qsample, ssample)), dist

    @cached_property
    def similarities(self):
        sims = dict(self._similarity_helper())
        self._samples = sorted(list(set(s for p in sims for s in p)))
        return sims | {frozenset([s]): 1 for s in self._samples}

    def get_dissimilarities(self):
        return {p: 1 - s for (p, s) in self.similarities.items()}

    def get_similarities(self):
        return self.similarities

    @property
    def samples(self):
        if self._samples is None:
            self.similarities
        return self._samples

    def _pair_dict_to_matrix(self, d):
        return np.vstack(
            [
                [
                    float(
                        d[frozenset([a, b])]
                    ) for b in self.samples
                ] for a in self.samples
            ]
        )

    def get_similarity_matrix(self):
        return self._pair_dict_to_matrix(self.get_similarities())

    def get_dissimilarity_matrix(self):
        return self._pair_dict_to_matrix(self.get_dissimilarities())

    @classmethod
    def from_filenames(cls, graph_fn, comparison_fns, *args, **kwargs):
        with open(graph_fn, "rb") as f:
            graph = pickle.load(f)
        return cls(graph, cls.load_pickles(comparison_fns), *args, **kwargs)
    
        
def main():
    args = handle_arguments()
    sim = SampleSimilarity.from_filenames(
        args.graph,
        tqdm(args.comparisons),
        args.samples
    )
    if args.out_type == "sim":
        print_mat(sim.get_similarity_matrix())
    elif args.out_type == "dis":
        print_mat(sim.get_dissimilarity_matrix())
    if args.embed:
        embed()

        
        
    

if __name__ == "__main__":
    main()
