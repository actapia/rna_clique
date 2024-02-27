#!/usr/bin/env bash
# This script further simplifies execution of the first phase of RNA-clique by
# automatically setting values for certain parameters.
#
# See do_filtering_step.sh to see further details about parameters for the first
# phase of RNA-clique.

function modify_help {
    sn=false
    while IFS= read -u 5 -r line; do
	if [[ -v last ]]; then
	    if [[ "$last" =~ arguments:$ ]] && [[ "$line" =~ ^$ ]]; then
		sn=true
	    else
		if [ "$sn" = false ]; then
		    echo "$last"
		fi
		sn=false
	    fi
	fi
	last="$line"
	#echo "$line hey"
    done 5< <(sed -ne '/arguments:$/,$ p' |
		  sed 's/\(^.*arguments\)/additional \1/' | grep -E -v "$1")
    if ! [[ "$last" =~ ^$ ]]; then
	echo "$last"
    fi
}

function modify_usage {
    read -r -a words
    re="$1"
    i=0
    #echo "$all_re"
    while [ "$i" -lt "${#words[@]}" ]; do
	if grep -E -q -e "$re" <<< "${words[$i]}"; then
	    #echo -n "${words[$i]} "
	    j=$((i+1))
	    if [[ "${words[$j]}" =~ ^[\]A-Z0-9\[]*$ ]] && \
		   ! grep -E -q -e "$3" <<< "${words[$j]}"; then
		#echo -n "${words[$j]}"
		i=$((i+1))
	    fi
	    #echo
	elif [ "$i" -ge "$2" ]; then
	    echo -n "${words[$i]} "
	fi
	i=$((i+1))
    done
    echo
}


declare -a rem
# n=50000
jobs="$(($(nproc) - 1))"
intermed=false
raw=false
do_help=false
readonly out_dir_flag="--output-dir"
readonly jobs_flag="--jobs"
readonly top_genes_flag="--top-genes"
readonly intermed_flag="--intermed"
readonly raw_flag="--raw"
readonly help_flag="--help"
declare -A ARG_SHORT
ARG_SHORT["$out_dir_flag"]="-o"
ARG_SHORT["$jobs_flag"]="-j"
ARG_SHORT["$top_genes_flag"]="-n"
ARG_SHORT["$intermed_flag"]="-i"
ARG_SHORT["$raw_flag"]="-r"
ARG_SHORT["$help_flag"]="-h"
declare -A ARG_HELP
ARG_HELP["$out_dir_flag"]="Parent output directory"
ARG_HELP["$jobs_flag"]="Number of parallel jobs to use"
ARG_HELP["$top_genes_flag"]="Number of top genes to select by k-mer coverage."
#ARG_HELP["$raw_flag"]="Keep the raw BLAST alignments."
#ARG_HELP["$intermed_flag"]="Keep intermediate results."
ARG_HELP["$help_flag"]="Print this help message and exit."
declare -A METAVAR
for flag in "$out_dir_flag" "$jobs_flag" "$top_genes_flag"; do
    mv="${flag//[^A-Za-z]/}"   
    METAVAR["$flag"]="${mv^^}"
done
declare -A REQUIRED
REQUIRED["$out_dir_flag"]=
REQUIRED["$top_genes_flag"]=
while [ "$#" -gt 0 ]; do
    case "$1" in
	"$out_dir_flag" | "${ARG_SHORT[$out_dir_flag]}")
	    shift;
	    out_dir="$1"
	    ;;
	"$jobs_flag" | "${ARG_SHORT[$jobs_flag]}")
	    shift
	    jobs="$1"
	    ;;
	"$top_genes_flag" | "${ARG_SHORT[$top_genes_flag]}")
	    shift
	    n="$1"
	    ;;
	"$intermed_flag" | "${ARG_SHORT[$intermed_flag]}")
	    intermed=true
	    ;;
	"$raw_flag" | "${ARG_SHORT[$raw_flag]}")
	    raw=true
	    ;;
	"$help_flag" | "${ARG_SHORT[$help_flag]}")
	    do_help=true
	    ;;
	*)
	    rem+=("$1")
	    ;;
    esac
    shift
done
declare -A EXCLUDE_ARGS
EXCLUDE_ARGS["--out-dir-1"]=
EXCLUDE_ARGS["--out-dir-2"]=
EXCLUDE_ARGS["--output-graph"]=
EXCLUDE_ARGS["--cache-dir"]=
EXCLUDE_ARGS["--keep-all"]=
EXCLUDE_ARGS["--help"]=
EXCLUDE_ARGS["--jobs"]=
EXCLUDE_ARGS["-n"]=
exclude_regex="$(printf "|%s" "${!EXCLUDE_ARGS[@]}")"
exclude_regex="(${exclude_regex:1})"
if [ "$do_help" = true ]; then
    # Print help.
    printf "Usage: $0 "
    longest=0
    for arg in "${!ARG_HELP[@]}"; do
	if [[ -v "ARG_SHORT[$arg]" ]]; then
	    arg_str="${ARG_SHORT[$arg]}"
	else	    
	    arg_str="$arg"
	fi
	if [[ -v "METAVAR[$arg]" ]]; then
	    arg_str="$arg_str ${METAVAR[$arg]}"
	fi
	if ! [[ -v "REQUIRED[$arg]" ]]; then
	    arg_str="[$arg_str]"
	fi
	printf "%s " "$arg_str"
	len="${#arg}"
	if [[ -v "ARG_SHORT[$arg]" ]]; then
	    len=$((len + "${#ARG_SHORT[$arg]}" + 2))
	fi
	if [ "$len" -gt "$longest" ]; then
	    longest="$len"
	fi
    done
    bash do_filtering_step.sh --help |
	grep '^Usage:' |
	modify_usage "$exclude_regex" 2 "DIR1"
    #printf "DIR [ ..."
    longest=$((longest+3))
    echo
    echo "positional arguments:"
    printf "  %-${longest}s" "DIR1 DIR2 ..."
    echo "Directories containing transcript FASTA files."
    echo
    echo "required arguments:"
    for arg in "${!ARG_HELP[@]}"; do
	if [[ -v "REQUIRED[$arg]" ]]; then
	    arg_str="$arg"
	    if [[ -v "ARG_SHORT[$arg]" ]]; then
		arg_str="${ARG_SHORT[$arg]}, $arg_str"
	    fi
	    printf "  %-${longest}s" "$arg_str"
	    echo "${ARG_HELP[$arg]}"
	fi
    done
    echo
    echo "main optional arguments:"
    for arg in "${!ARG_HELP[@]}"; do
	if ! [[ -v "REQUIRED[$arg]" ]]; then
	    arg_str="$arg"
	    if [[ -v "ARG_SHORT[$arg]" ]]; then
		arg_str="${ARG_SHORT[$arg]}, $arg_str"
	    fi
	    printf "  %-${longest}s" "$arg_str"
	    echo "${ARG_HELP[$arg]}"
	fi
    done

    echo
    #echo "additional arguments:"
    #echo "$exclude_regex"
    bash do_filtering_step.sh --help |
	sed -ne '/optional arguments:$/,$ p' |
	modify_help "$exclude_regex"
    exit 0
fi


if ! [[ -v out_dir ]]; then
    >&2 echo "Missing required argument --out-dir."
    exit 1
fi
if ! [[ -v n ]]; then
    >&2 echo "Missing required argument -n."
    exit 1
fi
intermed_add=()
if [ "$intermed" = true ]; then
    intermed_add=(--intermed-dir "$out_dir/intermed")
fi
raw_add=()
if [ "$raw" = true ]; then
    raw_add=(--raw-dir "$out_dir/raw")
fi
mkdir -p "$out_dir"
set -e
set -x
bash do_filtering_step.sh --jobs "$jobs" --out-dir-1 "$out_dir/od1" \
     --out-dir-2 "$out_dir/od2" -n "$n" --output-graph "$out_dir/graph.pkl" \
     --cache-dir "$out_dir/db_cache" "${intermed_add[@]}" "${raw_add[@]}" \
     --keep-all "${rem[@]}" 
