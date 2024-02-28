import argparse
import functools
import operator
import pickle
import re
from pathlib import Path
from find_homologs import eprint

from collections.abc import Iterable
from typing import Iterator

import pandas as pd
import networkx as nx

from tqdm import tqdm

# Since sum doesn't work on all objects.
sum_ = functools.partial(functools.reduce, operator.add)

def handle_arguments():
    parser = argparse.ArgumentParser(
        description="build a gene matches graph from gene matches tables"
    )
    parser.add_argument(
        "-i",
        "--inputs",
        nargs="+",
        type=Path,
        required=True,
        help="pickles containing the gene matches tables"
    )
    parser.add_argument(
        "-o",
        "--output-graph",
        type=Path,
        required=True,
        help="the output pickle of the gene matches graph"
    )
    return parser.parse_args()

def component_subgraphs(g : nx.Graph) -> Iterator[nx.Graph]:
    """Yields the connected components of the given graph as subgraphs."""
    for c in nx.connected_components(g):
        yield g.subgraph(c)

def make_edge(r):
    return (r[0], r[1]), (r[2], r[3])

def build_graph(dfs : Iterable[pd.DataFrame]) -> nx.Graph:
    """Build a gene matches graph from gene matches tables (dataframes).

    Each vertex in the gene matches graph is a tuple (s, g), where s is an
    identifier for the sample, and g is an identifier for a gene.

    There is an edge between two vertices (s1, g1) and (s2, g2) in the gene
    matches graph exactly when (s1, g1) and (s2, g2) appear in a gene matches
    table together with s1sample = s1, s1gene = g1, s2sample = s2, g2sample =
    g2, OR s1sample = s2, s1gene = g2, s2sample = s1, s2gene = g1.

    Parameters:
        dfs: The gene matches tables for the samples under consideration.

    Returns:
        The gene matches graph constructed from the given gene matches tables.
    """
    graph = nx.Graph()
    eprint("Building graph.")
    for df in dfs:
        all_cols = [[t + c for c in ["sample", "gene"]] for t in ["s", "q"]]
        for cols in all_cols:
            graph.add_nodes_from(
                tuple(s) for s in df[cols].itertuples(index=False)
            )
        graph.add_edges_from(
            make_edge(r) for r in df[sum_(all_cols)].itertuples(index=False)
        )
    return graph

def main():
    args = handle_arguments()
    graph = build_graph(pd.read_pickle(f) for f in tqdm(args.inputs))
    with open(args.output_graph, "wb") as f:
        pickle.dump(graph, f, pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    main()
