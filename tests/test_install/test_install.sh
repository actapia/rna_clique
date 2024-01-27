#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
readonly config_file="$DIR/minimal_config.yaml"
set -e
set -x
cd "$DIR"
cd ../../distance_sequence_simulator/
out_dir="$(mktemp -d)"
sim_data_dir="$out_dir/sim_data"
distance_dir="$out_dir/distance"
mkdir "$sim_data_dir"
mkdir "$distance_dir"

python generate_simulated_data.py -c "$config_file"  -j 0 -O "$sim_data_dir"
cd ..
bash typical_filtering_step.sh -n "$(python "$DIR/get_yaml_key.py" "$config_file" count)" -o "$distance_dir" "$sim_data_dir/"T*
PYTHONPATH="$(realpath "$DIR/../.."):$PYTHONPATH" python "$DIR/verify_distances.py" "$distance_dir" "$sim_data_dir/phylogeny.tree"
rm -r "$out_dir"
