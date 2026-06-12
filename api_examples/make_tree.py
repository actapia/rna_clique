import Bio.Phylo
import pandas as pd
import matplotlib as mpl
import nice_colorsys.rgb255

from pathlib import Path

from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from matplotlib import pyplot as plt
from nice_colorsys import rgb

from rna_clique import config as config_module
from rna_clique.viz.phylo_utils import (
    tril_jagged,
    draw_tree,
    get_clades,
    draw_clade_labels
)
from rna_clique.config import OutFileRule
from rna_clique.path_to_sample import path_to_sample as pts


def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager(
        description="Draw a phylogenetic tree for a given distance matrix."
    )
    arg_config.expose_fields_with_default_aliases(
        "matrix",
        required=True
    )
    arg_config.add_argument(
        "--metadata",
        "-M",
        type=Path,
        required=True,
        help="Path to metadata as comma-separated values."
    )
    arg_config.add_argument(
        "--group-by",
        "-g",
        required=True,
        nargs="+",
        help="Columns by which to group samples."
    )
    arg_config.add_argument(
        "--name-column",
        "-n",
        default="accession",
        help="Column providing name of sample."
    )
    arg_config.add_argument(
        "--output-plot",
        "-p",
        default={
            ("output_dir",): OutFileRule("output_dir", "tree.svg")
        }
    )
    arg_config.add_argument(
        "--output-tree",
        "-t",
        default={
            ("output_dir",): OutFileRule("output_dir", "nj_tree.tree")
        }
    )    
    arg_config.expose_fields_with_default_aliases("output_dir")
    return arg_config

def main():
    _, args, config = build_parser().get_arguments_and_config()
    sample_metadata = pd.read_csv(args.metadata)
    try:
        path_to_sample = {
            str(k): v for (k, v) in config.path_to_sample.items()
        }.__getitem__
    except AttributeError:
        path_to_sample = pts
    dis_df = pd.read_hdf(config.matrix).rename(
        index=path_to_sample,
        columns=path_to_sample,
    )
    nj_tree = DistanceTreeConstructor().nj(
        DistanceMatrix(
            list(dis_df.columns),
            tril_jagged(dis_df)
        )
    )
    nj_tree.root_at_midpoint()
    for c in nj_tree.get_nonterminals():
        c.name = None
    Bio.Phylo.write(nj_tree, args.output_tree, "newick")
    clades = dict(
        get_clades(
            nj_tree,
            sample_metadata,
            args.name_column,
            "genotype"
        )
    )
    clade_colors = {
        l: "#" + rgb(*x).to_rgb255().as_hex()
        for (l, x) in zip(clades, mpl.colormaps.get_cmap("tab10").colors)
    }
    draw_tree(nj_tree, clades=clades, colors=clade_colors)
    draw_clade_labels(plt.gca(), clades, clade_colors)
    plt.savefig(args.output_plot)

if __name__ == "__main__":
    main()
