import os
from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt

from rna_clique.viz.heatmap import draw_heatmap
from rna_clique.config import RNACliqueConfig

tutorial_doc_dir = Path(os.environ["RNA_CLIQUE"]) / "docs/tutorials/reads2tree"
rna_clique_out_dir = Path(os.environ["TUTORIAL_DIR"]) / "rna_clique_out"

def main():
    sample_metadata = pd.read_csv(tutorial_doc_dir / "tall_fescue_accs.csv")
    config = RNACliqueConfig.yaml_load(rna_clique_out_dir / "config.yaml")
    path_to_sample = {str(k): v for (k, v) in config.path_to_sample.items()}
    dis_df = pd.read_hdf(config.matrix).rename(
        index=path_to_sample.__getitem__,
        columns=path_to_sample.__getitem__,
    )
    draw_heatmap(
        dis_df,
        sample_metadata=sample_metadata,
        sample_name_column="accession",
        order_by="genotype",
        cmap="mako_r",
        digit_annot=2, # Show two digits of the distance.
        draw_group_labels=True, # Label according to genotype.
        label_padding_x = 0.05,
        label_padding_y = 0.05
    )
    plt.savefig(rna_clique_out_dir / "distance_heatmap.svg")
    
if __name__ == "__main__":
    main()

