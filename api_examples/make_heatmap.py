from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt

from rna_clique import config as config_module
from rna_clique.viz.heatmap import draw_heatmap
from rna_clique.config import OutFileRule
from rna_clique.path_to_sample import path_to_sample as pts

def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager(
        description="Draw a heatmap for a given distance matrix."
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
        help="Columns by which to group and order samples."
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
            ("output_dir",): OutFileRule("output_dir", "distance_heatmap.svg")
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
    draw_heatmap(
        dis_df,
        sample_metadata=sample_metadata,
        sample_name_column=args.name_column,
        order_by=args.group_by,
        cmap="mako_r",
        digit_annot=2, # Show two digits of the distance.
        draw_group_labels=True, # Label according to genotype.
        label_padding_x = 0.05,
        label_padding_y = 0.05
    )
    plt.savefig(args.output_plot)
    
if __name__ == "__main__":
    main()

