#!/usr/bin/env bash
set -e
set -x
declare -a rem
# n=50000
jobs="$(($(nproc) - 1))"
intermed=false
raw=false
while [ "$#" -gt 0 ]; do
    case "$1" in
	"--output-dir" | "-o")
	    shift;
	    out_dir="$1"
	    ;;
	"--jobs" | "-j")
	    shift
	    jobs="$1"
	    ;;
	"--top-genes" | "-n")
	    shift
	    n="$1"
	    ;;
	"--intermed" | "-i")
	    shift
	    intermed=true
	    ;;
	"--raw" | "-r")
	    shift
	    raw=true
	    ;;
	*)
	    rem+=("$1")
	    ;;
    esac
    shift
done
if ! [[ -v out_dir ]]; then
    >&2 echo "Missing required argument --out-dir."
    exit 1
fi
if ! [[ -v n ]]; then
    >&2 echo "Missing required argument -n."
    exit 1
fi
intermed_flag=()
if [ "$intermed" = true ]; then
    intermed_flag=(--intermed-dir "$out_dir/intermed")
fi
raw_flag=()
if [ "$raw" = true ]; then
    raw_flag=(--raw-dir "$out_dir/raw")
fi
mkdir -p "$out_dir"
set -e
set -x
bash do_filtering_step.sh --jobs "$jobs" --out-dir-1 "$out_dir/od1" --out-dir-2 "$out_dir/od2" -n "$n" --output-graph "$out_dir/graph.pkl" --cache-dir "$out_dir/db_cache" "${intermed_flag[@]}" "${raw_flag[@]}" --keep-all "${rem[@]}" 
