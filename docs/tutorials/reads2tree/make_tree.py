import os
import Bio.Phylo
import pandas as pd
import matplotlib as mpl

from pathlib import Path
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from matplotlib import pyplot as plt
from filtered_distance import SampleSimilarity
from phylo_utils import tril_jagged, draw_tree, get_clades, draw_clade_labels
from path_to_sample import path_to_sample
from make_subset import get_table_files
from nice_colorsys import *
from nice_colorsys.rgb255 import rgb255

tutorial_doc_dir = Path(os.environ["RNA_CLIQUE"]) / "docs/tutorials/reads2tree"
rna_clique_out_dir = Path(os.environ["TUTORIAL_DIR"]) / "rna_clique_out"

def main():
    sample_metadata = pd.read_csv(tutorial_doc_dir / "tall_fescue_accs.csv")
    similarity_computer = SampleSimilarity.from_filenames(
        rna_clique_out_dir / "graph.pkl",
        list(get_table_files(rna_clique_out_dir / "od2"))
    )
    nj_tree = DistanceTreeConstructor().nj(
        DistanceMatrix(
            [path_to_sample(p) for p in similarity_computer.samples],
            tril_jagged(similarity_computer.get_dissimilarity_matrix())
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
