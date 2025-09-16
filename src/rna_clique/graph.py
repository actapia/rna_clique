import functools

import networkx as nx

from typing import Iterator

def connected_component_subgraphs(
        g,
        connected_components
) -> Iterator:
    """Yields the connected components of the given graph as subgraphs."""
    for c in connected_components(g):
        yield g.subgraph(c)

component_subgraphs = functools.partial(
    connected_component_subgraphs,
    connected_components=nx.connected_components
)
weak_component_subgraphs = functools.partial(
    connected_component_subgraphs,
    connected_components=nx.weakly_connected_components
)
