import itertools
from pysat.examples.rc2 import RC2
from collections import defaultdict

import networkx as nx
from sympy.logic.boolalg import to_cnf, to_int_repr
from sympy.abc import A, B, C

from pysat.formula import WCNF

eq = ~(A ^ B)
neq = A ^ B
edge_eq = to_cnf((C >> eq) & (eq >> C), True)
edge_neq = to_cnf((C >> neq) & (neq >> C), True)

def id_dicts(n=1, start=0):
    def incr():
        nonlocal cnt
        cnt += 1
        return cnt
    cnt = start - 1
    return [defaultdict(incr) for _ in range(n)]

def cnf_subs(expr, a, b, edge):
    replace = {1: a, 2: b, 3: edge}
    replace |= {-k: -v for (k, v) in replace.items()}
    return [[replace[k] for k in c] for c in to_int_repr(expr.args, [A, B, C])]

def to_maxsat_problem(g):
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

def sign(x):
    return x // abs(x)

def apply_assignment(g, assignment, dec_nodes):
    for node, id_ in dec_nodes.items():
        g.nodes[node]["strand"] = sign(assignment[id_ - 1])

def sat_assign_strands(g):
    wcnf, dec_nodes, _ = to_maxsat_problem(g)
    with RC2(wcnf) as rc2:
        rc2.compute()
        cost = rc2.cost
        apply_assignment(g, next(iter(rc2.enumerate())), dec_nodes)
        return cost

def example_graph():
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
