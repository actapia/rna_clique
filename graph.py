import functools

import networkx as nx

from typing import Iterator, TypeVar, Callable

T = TypeVar("T")

def connected_component_subgraphs(
        g: T,
        connected_components: Callable[[T], set]
) -> Iterator[T]:
    """Returns iterator of connected components of the given graph as subgraphs.

    Parameters:
        g:                    The NetworkX Graph-like object.
        connected_components: Function to get connected components from g.
    """
    return map(g.subgraph, connected_components(g))

component_subgraphs = functools.partial(
    connected_component_subgraphs,
    connected_components=nx.connected_components
)
weak_component_subgraphs = functools.partial(
    connected_component_subgraphs,
    connected_components=nx.weakly_connected_components
)
