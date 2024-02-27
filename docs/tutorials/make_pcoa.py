import os
from pathlib import Path

import skbio as skb
import pandas as pd
from matplotlib import pyplot as plt

from IPython import embed
from filtered_distance import SampleSimilarity
from path_to_sample import path_to_sample

tutorial_doc_dir = Path(os.environ["RNA_CLIQUE"]) / "docs" / "tutorials"
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
    embed()
    # 3D PCoA
    pcoa_results = skb.stats.ordination.pcoa(
        skb.DistanceMatrix(dis_df, ids=dis_df.columns)
    )
    pcoa_results.plot(
        df=sample_metadata.set_index("accession"),
        column="genotype",
    )
    plt.savefig(rna_clique_out_dir / "pcoa_3d.svg")
    # 2D PCoA
    pcoa_results_2d = skb.stats.ordination.pcoa(
        skb.DistanceMatrix(dis_df, ids=dis_df.columns),
        number_of_dimensions=2
    )
    plt.figure()
    for g, df in sample_metadata.join(
            pcoa_results_2d.samples[["PC1","PC2"]],
            "accession"
    ).groupby("genotype"):
        plt.scatter(df["PC1"], df["PC2"], label=g)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    plt.savefig(rna_clique_out_dir / "pcoa_2d.svg")
    

if __name__ == "__main__":
    main()

