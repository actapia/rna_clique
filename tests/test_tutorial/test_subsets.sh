#!/usr/bin/env bash
set -e
cd rna_clique
export RNA_CLIQUE="$PWD"
cd ..
[ -d "tutorial" ]
cd tutorial
export TUTORIAL_DIR="$PWD"
cd "$RNA_CLIQUE"
eval "$("$HOME/miniconda3/bin/conda" shell.bash hook)"
conda activate rna-clique
python make_subset.py -I "$TUTORIAL_DIR"/rna_clique_out/config.yaml \
                      -O "$TUTORIAL_DIR"/infected_subset_out \
                      -y SRR2321388 SRR8003761 SRR7990321 SRR8003736
! [ -f "tutorial/infected_subset_out/distance_matrix.h5" ]
python filtered_distance.py -O "$TUTORIAL_DIR"/infected_subset_out
[ -f "tutorial/infected_subset_out/distance_matrix.h5" ]
components="$(python plot_component_sizes.py --statistics \
                     -A "$TUTORIAL_DIR"/infected_subset_out)"
[ "$components" -ge 10000 ]

