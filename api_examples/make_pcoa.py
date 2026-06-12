import os
import argparse
import functools

import pandas as pd
import matplotlib as mpl

from pathlib import Path

from matplotlib import pyplot as plt

from rna_clique import config as config_module
from rna_clique.viz import pcoa
from rna_clique.config import OutFileRule
from rna_clique.path_to_sample import path_to_sample as pts

def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager(
        description="Draw a PCoA plot for a given distance matrix."
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
            ("output_dir",): OutFileRule("output_dir", "pcoa.svg")
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
    draw_pcoa = functools.partial(
        pcoa.draw_pcoa,
        dis_df,
        sample_metadata,
        group_by=args.group_by,
        sample_name_column=args.name_column,
        colors=mpl.colormaps.get_cmap("tab10"),
        edgecolors="black",
        linewidth=0.3,
    )
    # 2D PCoA
    draw_pcoa(dimensions=2)
    plt.savefig(args.output_plot)

if __name__ == "__main__":
    main()

