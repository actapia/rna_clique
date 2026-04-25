import os

import Bio.Phylo
import pandas as pd
import matplotlib as mpl
import nice_colorsys.rgb255

from pathlib import Path

from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from matplotlib import pyplot as plt
from nice_colorsys import rgb

from phylo_utils import tril_jagged, draw_tree, get_clades, draw_clade_labels
from config import RNACliqueConfig

tutorial_doc_dir = Path(os.environ["RNA_CLIQUE"]) / "docs/tutorials/reads2tree"
rna_clique_out_dir = Path(os.environ["TUTORIAL_DIR"]) / "rna_clique_out"

def main():
    sample_metadata = pd.read_csv(tutorial_doc_dir / "tall_fescue_accs.csv")
    config = RNACliqueConfig.yaml_load(rna_clique_out_dir / "config.yaml")
    nj_tree = DistanceTreeConstructor().nj(
        DistanceMatrix(
            list(config.path_to_sample.values()),
            tril_jagged(pd.read_hdf(config.matrix))
        )
    )
    nj_tree.root_at_midpoint()
    for c in nj_tree.get_nonterminals():
        c.name = None
    Bio.Phylo.write(nj_tree, rna_clique_out_dir / "nj_tree.tree", "newick")
    clades = dict(get_clades(nj_tree, sample_metadata, "accession", "genotype"))
    clade_colors = {
        l: "#" + rgb(*x).to_rgb255().as_hex()
        for (l, x) in zip(clades, mpl.colormaps.get_cmap("tab10").colors)
    }
    draw_tree(nj_tree, clades=clades, colors=clade_colors)
    draw_clade_labels(plt.gca(), clades, clade_colors)
    plt.savefig(rna_clique_out_dir / "nj_tree.svg")

if __name__ == "__main__":
    main()
