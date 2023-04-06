#!/usr/bin/env bash
if ! [[ -v INCSEQ ]]; then
    INCSEQ="./incseq.sh"
fi
jobs=$(($(nproc) - 1))
step=1000
transcripts="transcripts.fasta"
clean_level=0
OUTDIR_ROOT="."
declare -a fasta_files
while [ "$#" -gt 0 ]; do
    case "$1" in
	"-j" | "--jobs")
	    shift
	    jobs="$1"
	    ;;
	"-m" | "--min")
	    shift
	    min="$1"
	    ;;
	"-s" | "--step")
	    shift
	    step="$1"
	    ;;
	"-M" | "--max")
	    shift
	    max="$1"
	    ;;
	"-t" | "--transcripts")
	    shift
	    transcripts="$1"
	    ;;
	"-c" | "--clean-level")
	    shift
	    clean_level="$1"
	    ;;
	"-O" | "--out")
	    shift
	    OUTDIR_ROOT="$1"
	    ;;
	*)
	    fasta_files+=("$1")
	    ;;
    esac
    shift
done
set -e
set -x
if ! [[ -v min ]]; then
    min="$step"
fi
if ! [[ -v max ]]; then
    max=0
    for f in "${fasta_files[@]}"; do
	gene_count="$(perl select_top_genes/count_genes.pl -t "$f/$transcripts")"
	if [ "$gene_count" -gt "$max" ]; then
	    max="$gene_count"
	fi
    done
fi
export clean_level
export OUTDIR_ROOT
>&2 echo "Incseq: $INCSEQ"
>&2 echo "Running $("$INCSEQ" "$min" "$step" "$max" | wc -l) tests with $jobs jobs."
"$INCSEQ" "$min" "$step" "$max" | TRY_PARAMS_JOBS="$jobs" bash try_params.sh -n "${fasta_files[@]}" --keep-all
