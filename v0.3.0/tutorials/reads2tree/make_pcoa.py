import os
import argparse
import functools

import pandas as pd
import matplotlib as mpl

from pathlib import Path

from matplotlib import pyplot as plt

from rna_clique.viz import pcoa
from rna_clique.config import RNACliqueConfig

tutorial_doc_dir = Path(os.environ["RNA_CLIQUE"]) / "docs/tutorials/reads2tree"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("rna_clique_out", type=Path)
    args = parser.parse_args()
    sample_metadata = pd.read_csv(tutorial_doc_dir / "tall_fescue_accs.csv")
    config = RNACliqueConfig.yaml_load(args.rna_clique_out / "config.yaml")
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
    plt.savefig(args.rna_clique_out / "pcoa_2d.svg")    
    # 3D PCoA
    draw_pcoa(dimensions=3)
    plt.savefig(args.rna_clique_out / "pcoa_3d.svg")

if __name__ == "__main__":
    main()

