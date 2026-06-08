#!/usr/bin/env bash
set -e
cd rna_clique
export RNA_CLIQUE="$PWD"
cd ..
[ -d "tutorial" ]
cd tutorial
export TUTORIAL_DIR="$PWD"
cd "$RNA_CLIQUE"
. rna_clique_venv/bin/activate
python -m rna_clique.make_subset -I "$TUTORIAL_DIR"/rna_clique_out/config.yaml \
                                 -O "$TUTORIAL_DIR"/infected_subset_out \
                                 -y SRR2321388 SRR8003761 SRR7990321 SRR8003736
! [ -f "$TUTORIAL_DIR/infected_subset_out/distance_matrix.h5" ]
python -m rna_clique.filtered_distance -O "$TUTORIAL_DIR"/infected_subset_out
[ -f "$TUTORIAL_DIR/infected_subset_out/distance_matrix.h5" ]
components="$(python rna_clique.plot_component_sizes --statistics m \
                     -A "$TUTORIAL_DIR"/infected_subset_out | cut -d' ' -f4)"
[ "$components" -ge 10000 ]

