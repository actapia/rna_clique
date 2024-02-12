#!/usr/bin/env bash
{ read -r a; read -r b; } < <(bash do_filtering_step.sh "$@" | tail -n2)
python plot_component_sizes.py --statistics=m "$a" --samples "$b"
