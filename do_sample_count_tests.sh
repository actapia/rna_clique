#!/usr/bin/env bash

readonly MIN_TEST=4

# get_seeded_random function from
# https://www.gnu.org/software/coreutils/manual/html_node/Random-sources.html#Random-sources
#
# Thanks to Charles Duffy on Stack Overflow for finding this.
# https://stackoverflow.com/a/41962458
function get_seeded_random
{
    seed="$1"
    openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
	    </dev/zero 2>/dev/null
}

if ! [[ -v INCSEQ ]]; then
    INCSEQ="./incseq.sh"
fi
declare -a fasta_files
declare -a other
out_root="."
out_fn="top_n_tests_%d_samples"
do_shuffle=true
while [ "$#" -gt 0 ]; do
    case "$1" in
	"-O" | "--out-dir")
	    shift
	    out_root="$1"
	    ;;
	"-o" | "--out")
	    shift
	    out_fn="$1"
	    ;;
	"--no-shuffle")
	    do_shuffle=false;
	    ;;
	-*)
	    other+=("$1")
	    shift
	    other+=("$1")
	    ;;
	*)
	    fasta_files+=("$1")
	    ;;
    esac
    shift
done
# set -e
# set -x
declare -a shuffled
export -f get_seeded_random
if [ "$do_shuffle" = true ]; then
    readarray -t shuffled < <(shuf --random-source=<(get_seeded_random 485) -e "${fasta_files[@]}")
else
    shuffled=("${fasta_files[@]}")
fi
declare -a current
# for f in "${shuffled[@]}"; do
#     echo "$f"
# done
for f in "${shuffled[@]}"; do
    current+=("$f")
    if [ "${#current[@]}" -ge 4 ]; then
	out="$out_root/out_${#current[@]}_samples"
	results="$out_root/$(printf "$out_fn" "${#current[@]}")"
	#echo "${other[@]}" -O "$out" "${current[@]}" "|" tee "$results"
	bash do_top_n_tests.sh "${other[@]}" -O "$out" "${current[@]}" | tee "$results"
    fi
done
