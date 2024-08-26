import os
from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt

from filtered_distance import SampleSimilarity
from path_to_sample import path_to_sample
from heatmap import draw_heatmap

tutorial_doc_dir = Path(os.environ["RNA_CLIQUE"]) / "docs/tutorials/reads2tree"
rna_clique_out_dir = Path(os.environ["TUTORIAL_DIR"]) / "rna_clique_out"

def main():
    sample_metadata = pd.read_csv(tutorial_doc_dir / "tall_fescue_accs.csv")
    similarity_computer = SampleSimilarity.from_filenames(
        rna_clique_out_dir / "graph.pkl",
        (rna_clique_out_dir / "od2").glob("*.pkl")
    )
    dis_df = similarity_computer.get_dissimilarity_df().rename(
        index=path_to_sample,
        columns=path_to_sample,
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

