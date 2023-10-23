import Bio.Phylo
import numpy as np

def tril_jagged(mat):
    return [r[:(i+1)].tolist() for (i, r) in enumerate(np.tril(mat))]

def draw_tree(tree):
    ax = plt.axes()
    Bio.Phylo.draw(tree, do_show=False, axes=ax)
    ax.axis("off")