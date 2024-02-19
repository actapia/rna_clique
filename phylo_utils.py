import Bio.Phylo
import numpy as np
from numbers import Number
from matplotlib import pyplot as plt
from Bio.Phylo import BaseTree

def tril_jagged(mat: np.ndarray) -> list[list[Number]]:
    """Get lower triangle of a matrix as a jagged 2-dimensional list."""
    return [r[:(i+1)].tolist() for (i, r) in enumerate(np.tril(mat))]

def draw_tree(tree: BaseTree):
    ax = plt.axes()
    Bio.Phylo.draw(tree, do_show=False, axes=ax)
    ax.axis("off")
