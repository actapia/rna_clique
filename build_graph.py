import functools
import operator
import pickle
import config as config_module
import app
from app import eprint, set_except_hook
from gene_matches_tables import get_table_files, read_table

from collections.abc import Iterable
from pathlib import Path

import pandas as pd
import networkx as nx

from tqdm import tqdm

# Since sum doesn't work on all objects.
sum_ = functools.partial(functools.reduce, operator.add)

def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager(
        description=(
            "Build the gene matches graph from the gene matches tables."
        ),
    )
    arg_config.expose_fields_with_default_aliases(
        "tables_dir",
        "graph",
        required=True,
    )
    arg_config.expose_fields_with_default_aliases(
        "output_dir",
    )
    arg_config.add_output_config_argument()
    return arg_config

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
    with set_except_hook():
        _, args, config = build_parser().get_arguments_and_config()
    with set_except_hook(config.verbose):
        config_module.RNACliqueConfig.validate_dir(config.tables_dir)
        tables = list(get_table_files(config.tables_dir))
        if not tables:
            eprint(
                "Warning: No gene matches tables found in {}".format(
                    config.tables_dir
                )
            )
        graph = build_graph(
            read_table(f) for f in tqdm(tables)
        )
        with open(config.graph, "wb") as f:
            pickle.dump(graph, f, pickle.HIGHEST_PROTOCOL)
        config.mark_finish()
        if args.output_config:
            config.yaml_save(args.output_config)

if __name__ == "__main__":
    main()
