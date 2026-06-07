import itertools

import networkx as nx
import sympy.logic.boolalg

from collections import defaultdict
from collections.abc import Mapping
from typing import Any

from pysat.examples.rc2 import RC2
from sympy.logic.boolalg import to_cnf, to_int_repr
from sympy.abc import A, B, C
from pysat.formula import WCNF

eq = ~(A ^ B)
neq = A ^ B
edge_eq = to_cnf((C >> eq) & (eq >> C), True)
edge_neq = to_cnf((C >> neq) & (neq >> C), True)

def id_dicts(n: int = 1, start: int = 0) -> list[defaultdict[Any, int]]:
    """Return multiple defaultdicts using a shared counter value as defaults.

    The returned defaultdicts have the property that no two items in any of the
    dicts share the same value. This effectively assigns an ID to every dict
    item that is unique among the created defaultdicts. Note that this means
    that two items with the same key will never share an ID. For example, if we
    created two defaultdicts, a and b, a["foo"] != b["foo"] because the items
    are in different dicts.

    Parameters:
        n (int):     Number of defaultdicts with shared counter to create.
        start (int): Initial value for the counter.

    Returns:
        n defaultdicts whose default values are provided by a shared counter.
    """
    def incr() -> int:
        nonlocal cnt
        cnt += 1
        return cnt
    cnt = start - 1
    return [defaultdict(incr) for _ in range(n)]

# TODO: This could be made more generic, and possibly simpler, by accepting an
# arbitrary number of arguments.
def cnf_subs(
        expr: sympy.logic.boolalg.And,
        a: int,
        b: int,
        edge: int,
) -> list[list[int]]:
    """Get integer representation of CNF formula with variables substituted.

    This function accepts a CNF formula and the variables to substitute,
    represented by integers. The CNF formula should have three variables,
    representing the two endpoints of an edge and the edge itself.  The returned
    value is a list of lists representing a conjunction of disjunctions of
    variables. Negated variables are represented by their additive inverses.

    Parameters:
        expr:       CNF formula with endpoints and edge as variables.
        a (int):    Variable representing one endpoint of the edge.
        b (int):    Variable representing the other endpoint of the edge.
        edge (int): Variable representing an edge.

    Returns:
        A integer representation of expr, with a, b, and edge substituted in.
    """
    replace = {1: a, 2: b, 3: edge}
    replace |= {-k: -v for (k, v) in replace.items()}
    return [[replace[k] for k in c] for c in to_int_repr(expr.args, [A, B, C])]

def to_maxsat_problem(
        g: nx.Graph
) -> tuple[
    WCNF,
    Mapping[tuple[str, int, int], int],
    Mapping[tuple[str, int, int], int]
]:
    """Create a MaxSAT problem for assigning the given strand graph.

    This function returns a formula in weighted conjunctive normal form (WCNF)
    such that an optimal assignment of values to the variables appearing in the
    formula corresponds to an optimal assignment of strands to nodes in the
    strand graph.

    In this formulation of the problem in MaxSAT, each vertex and each edge
    corresponds to a MaxSAT variable. A vertex's variable is true if and only if
    the vertex should be assigned 1. An edge's variable is true if and only if
    the edge is satisfied.

    The hard clauses---those that must be satisfied by all solutions---enforce
    the requirement that an edge's variable is true if and only if the edge is
    satisfied. When the edge has weight -1, the corresponding clauses state that
    the edge's variable can be satisfied if and only if its endpoints are not
    equal. When the edge has weight 1, the corresponding clauses state that the
    edge's variable can be satisfied if and only if its edgepoints are equal.

    The soft clauses---those that do not have to be satisfied by every solution,
    but which increase the score of the solution when satisfied---are simply the
    variables corresponding to the edges.

    Thus, the MaxSAT solver tries to maximize the number of edge variables set
    to true without setting an edge variable to true when the edge is not
    actually satisfied.

    In addition to the WCNF formula, this function returns two Mappings. The
    first maps vertices to MaxSAT variables. The second maps edges to MaxSAT
    variables.

    Parameters:
        g: The strand graph for which to make a MaxSAT problem.

    Returns:
        The MaxSAT WNCF formula and mappings from nodes and edges to variables.
    """
  
    form = WCNF()
    node_to_var, edge_to_var = id_dicts(2, 1)
    for edge in g.edges:
        expr = edge_eq if g.edges[edge]["weight"] == 1 else edge_neq
        for clause in cnf_subs(
                expr,
                node_to_var[edge[0]],
                node_to_var[edge[1]],
                edge_to_var[edge]
        ):
            form.append(clause)
        form.append([edge_to_var[edge]], weight=1)
    return (
        form,
        node_to_var,
        edge_to_var,
    )

def sign(x: int) -> int:
    """Return the sign of the integer, -1 if x is negative or 1 if positive.

    Naturally, this function raises a ValueError if x is 0.

    Parameters:
        x (int): The nonzero integer for which to obtain a sign.

    Returns:
        The sign of x, -1 if x is negative, or 1 if x is positive.
    """
    return x // abs(x)

def apply_assignment(
        g: nx.Graph,
        assignment: Mapping[int, int],
        dec_nodes: Mapping[tuple[str, int, int], int]
):
    """Apply a MaxSAT variable assignment to a strand graph.

    This function accepts the strand graph to modify, an assignment of MaxSAT
    variables, and a Mapping from strand graph nodes to their corresponding
    MaxSAT variables.

    The assignment is expected to be a Mapping from variables (represented by
    integers) to their assignments (also represented by integers). If the value
    for a given variable key is negative, the variable is assumed to be assigned
    false in the MaxSAT solution; this corresponds to an assignment of -1 in the
    strand graph. Likewise, if the value for a given variable key is positive,
    the variable is assumed to be assigned true in the MaxSAT solution, and this
    is interpreted as an assignment of +1 to the strand graph node.

    Parameters:
        g:          The graph for which to assign strands to nodes.
        assignment: The MaxSAT solution assigning truth values to variables.
        dec_nodes:  Mapping from graph nodes to MaxSAT variables.
    """
    for node, id_ in dec_nodes.items():
        g.nodes[node]["strand"] = sign(assignment[id_ - 1])

def sat_assign_strands(g: nx.Graph) -> int:
    """Assign strands to nodes optimally using MaxSAT.

    Parameters:
        g: The graph for which to assign strands to nodes.

    Returns:
        The cost of the MaxSAT solution, the number of violated edges.
    """
    wcnf, dec_nodes, _ = to_maxsat_problem(g)
    with RC2(wcnf) as rc2:
        rc2.compute()
        cost = rc2.cost
        apply_assignment(g, next(iter(rc2.enumerate())), dec_nodes)
        return cost

def example_graph() -> nx.Graph:
    """Make a very simple example strand graph for testing strand assignment.

    The returned graph is the complete graph with three vertices. It is
    inconsistent because exactly one of the edges has weight -1.

    Parameters:
        An example graph for testing strand assignment.
    """
    g = nx.Graph()
    nodes = [f"sample_{l}" for l in "ABC"]
    edges = list(itertools.combinations(nodes, r=2))
    g.add_nodes_from(nodes)
    g.add_edges_from(edges)
    g.edges[edges[0]]["weight"] = 1
    g.edges[edges[1]]["weight"] = -1
    g.edges[edges[2]]["weight"] = 1
    return g

def main():
    g = example_graph()
    wcnf, dec_nodes, dec_edges = to_maxsat_problem(g)
    with RC2(wcnf) as rc2:
        rc2.compute()
        from IPython import embed; embed()

if __name__ == "__main__":
    main()
