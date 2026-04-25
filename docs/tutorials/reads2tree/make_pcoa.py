import os
import functools

import pandas as pd
import matplotlib as mpl

import pcoa

from pathlib import Path

from matplotlib import pyplot as plt

from config import RNACliqueConfig

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
    draw_pcoa = functools.partial(
        pcoa.draw_pcoa,
        dis_df,
        sample_metadata,
        group_by="genotype",
        sample_name_column="accession",
        colors=mpl.colormaps.get_cmap("tab10"),
        edgecolors="black",
        linewidth=0.3,
    )
    # 2D PCoA
    draw_pcoa(dimensions=2)
    plt.savefig(rna_clique_out_dir / "pcoa_2d.svg")    
    # 3D PCoA
    draw_pcoa(dimensions=3)
    plt.savefig(rna_clique_out_dir / "pcoa_3d.svg")

if __name__ == "__main__":
    main()

