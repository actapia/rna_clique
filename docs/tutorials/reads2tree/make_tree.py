import os
import Bio.Phylo

from pathlib import Path
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from matplotlib import pyplot as plt
from filtered_distance import SampleSimilarity
from phylo_utils import tril_jagged, draw_tree
from path_to_sample import path_to_sample
from make_subset import get_table_files

rna_clique_out_dir = Path(os.environ["TUTORIAL_DIR"]) / "rna_clique_out"

def main():
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
    draw_tree(nj_tree)
    plt.savefig(rna_clique_out_dir / "nj_tree.svg")
    
if __name__ == "__main__":
    main()
