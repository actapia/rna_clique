import argparse
import functools
import operator
import pickle
from pathlib import Path
from find_homologs import eprint

import pandas as pd
import networkx as nx

from tqdm import tqdm

from IPython import embed

sum_ = functools.partial(functools.reduce, operator.add)

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputs", nargs="+", type=Path, required=True)
    parser.add_argument("-o", "--output-graph", type=Path)
    return parser.parse_args()

def make_edge(r):
    return ((r[0], r[1]), (r[2], r[3]))

def build_graph(dfs):
    graph = nx.Graph()
    eprint("Building graph.")
    for df in dfs:
        all_cols = [[t + c for c in ["sample", "gene"]] for t in ["s", "q"]]
        for cols in all_cols:
            graph.add_nodes_from(tuple(s) for s in df[cols].itertuples(index=False))
        graph.add_edges_from(make_edge(r) for r in df[sum_(all_cols)].itertuples(index=False))
    return graph
        
def main():
    args = handle_arguments()
    graph = build_graph(pd.read_pickle(f) for f in tqdm(args.inputs))
    if args.output_graph:
        with open(args.output_graph, "wb") as f:
            pickle.dump(graph, f, pickle.HIGHEST_PROTOCOL)
    else:
        embed()

if __name__ == "__main__":
    main()
