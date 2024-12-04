#!/usr/bin/env bash
#
# This script performs the first "phase" of RNA-clique, which involves the
# following steps:
#
# 1. Selecting the top n genes from each sample.
# 2. Obtaining the gene matches tables for each pair of samples. (This involves
#    a pairwise reciprocal BLAST search.)
# 3. Building the gene matches graph from the gene matches tables.
#
# This script allows the user to set various parameters that control the
# behavior of RNA-clique in this first phase. For many parameters, reasonable
# defaults are provided, so it is not always necessary to provide every
# parameter.
readonly top_n_flag="-n"
readonly top_N_flag="-N"
readonly transcripts_flag="--transcripts"
readonly out_dir_1_flag="--out-dir-1"
readonly out_dir_2_flag="--out-dir-2"
readonly pattern_flag="--pattern"
readonly evalue_flag="--evalue"
readonly sample_regex_flag="--sample-regex"
readonly keep_all_flag="--keep-all"
readonly output_graph_flag="--output-graph"
readonly gene_regex_flag="--gene-regex"
readonly jobs_flag="--jobs"
readonly cache_dir_flag="--cache-dir"
readonly help_flag="--help"
# Default values.
transcripts_fn="transcripts.fasta"
top_n=10000
top_N=1
pattern="^.*cov_([0-9]+(?:\.[0-9]+))_g([0-9]+)_i([0-9]+)"
gene_regex="^.*g([0-9]+)_i([0-9]+)"
# SAMPLE_RE may also be provided as an environment variable.
if [[ -v SAMPLE_RE ]]; then
    sample_regex="$SAMPLE_RE"
else
    sample_regex="^(.*?)_.*$"
fi
do_help=false
evalue=1e-99
keep_all=false
jobs=1
declare -A ARG_HELP
ARG_HELP["$top_n_flag"]="Top n genes to select."
ARG_HELP["$top_N_flag"]="Count a match if it is within top N in both directions."
ARG_HELP["$transcripts_flag"]="Name of transcripts files."
ARG_HELP["$out_dir_1_flag"]="Intermediate output file directory containing top genes."
ARG_HELP["$out_dir_2_flag"]="Intermediate output file directory containing matches."
ARG_HELP["$pattern_flag"]="Regular expression for getting coverage values."
ARG_HELP["$evalue_flag"]="Cutoff evalue to use in BLAST searches."
ARG_HELP["$sample_regex_flag"]="Regular expression matching sample names."
ARG_HELP["$gene_regex_flag"]="Regular expression matching gene and isoform IDs."
ARG_HELP["$keep_all_flag"]="Keep all matches in case of a tie."
ARG_HELP["$output_graph_flag"]="Path to output graph."
ARG_HELP["$jobs_flag"]="Number of parallel jobs to use."
ARG_HELP["$cache_dir_flag"]="Directory in which to store BLAST DBs."
ARG_HELP["$help_flag"]="Show this help and exit."
declare -A METAVAR
for flag in "$transcripts_flag" "$out_dir_1_flag" "$out_dir_2_flag" \
	       "$pattern_flag" "$evalue_flag" "$sample_regex_flag" \
	       "$gene_regex_flag" "$output_graph_flag" "$jobs_flag" \
	       "$cache_dir_flag" "$help_flag"; do
    mv="${flag//[^A-Za-z]/}"   
    METAVAR["$flag"]="${mv^^}"
done
METAVAR["$top_n_flag"]="TOPNGENES"
METAVAR["$top_N_flag"]="TOPNMATCHES"
declare -A REQUIRED
for flag in 	"$out_dir_1_flag" \
		"$out_dir_2_flag" \
		"$output_graph_flag"; do
    REQUIRED["$flag"]=
done
declare -a dirlist
while [ "$#" -gt 0 ]; do
    case "$1" in
	"$top_n_flag")
	    shift;
	    top_n="$1"
	    ;;
	"$top_N_flag")
	    shift;
	    top_N="$1"
	    ;;
	"$transcripts_flag")
	    shift;
	    transcripts_fn="$1"
	    ;;
	"$out_dir_1_flag")
	    shift;
	    out_dir_1="$1"
	    ;;
	"$out_dir_2_flag")
	    shift;
	    out_dir_2="$1"
	    ;;
	"$pattern_flag")
	    shift;
	    pattern="$1"
	    ;;
	"$evalue_flag")
	    shift;
	    evalue="$1"
	    ;;
	"$sample_regex_flag")
	    shift;
	    sample_regex="$1"
	    ;;
	"$gene_regex_flag")
	    shift;
	    gene_regex="$1"
	    ;;
	"$keep_all_flag")
	    keep_all=true
	    ;;
	"$output_graph_flag")
	    shift;
	    output_graph="$1"
	    ;;
	"$jobs_flag")
	    shift;
	    jobs="$1"
	    ;;
	"$cache_dir_flag")
	    shift;
	    cache_dir="$1"
	    ;;
	"$help_flag")
	    do_help=true
	    ;;
	*)
	    if ! [ -d "$1" ]; then
		>&2 echo "Directory $1 does not exist."
		exit 1
	    fi
	    dirlist+=("$1")
	    ;;
    esac
    shift
done

if [ "$do_help" = true ]; then
    # Print help.
    printf "Usage: $0 "
    longest=0
    for arg in "${!ARG_HELP[@]}"; do
	arg_str="$arg"
	if [[ -v "METAVAR[$arg]" ]]; then
	    arg_str="$arg_str ${METAVAR[$arg]}"
	fi
	if ! [[ -v "REQUIRED[$arg]" ]]; then
	    arg_str="[$arg_str]"
	fi
	printf "%s " "$arg_str"
	if [ "${#arg}" -gt "$longest" ]; then
	    longest="${#arg}"
	fi
    done
    printf "DIR1 DIR2 ..."
    longest=$((longest+3))
    echo
    echo
    echo "positional arguments:"
    printf "  %-${longest}s" "DIR1 DIR2 ..."
    echo "Directories containing transcript FASTA files."
    echo
    echo "required arguments:"
    for arg in "${!ARG_HELP[@]}"; do
	if [[ -v "REQUIRED[$arg]" ]]; then
	    printf "  %-${longest}s" "$arg"
	    echo "${ARG_HELP[$arg]}"
	fi
    done
    echo
    echo "optional arguments:"
    for arg in "${!ARG_HELP[@]}"; do
	if ! [[ -v "REQUIRED[$arg]" ]]; then
	    printf "  %-${longest}s" "$arg"
	    echo "${ARG_HELP[$arg]}"
	fi
    done
    exit 0
fi

set -x

if ! [[ -v top_n ]] ; then
    >&2 echo "Missing required argument $top_n_flag."
    exit 1
fi
if ! [[ -v top_N ]]; then
    >&2 echo "Missing required argument $top_N_flag."
    exit 1
fi
if ! [[ -v transcripts_fn ]]; then
    >&2 echo "Missing required argument $transcripts_flag."
    exit 1
fi
if ! [[ -v out_dir_1 ]]; then
    >&2 echo "Missing required argument $out_dir_1_flag."
    exit 1
fi
if ! [[ -v out_dir_2 ]]; then
    >&2 echo "Missing required argument $out_dir_2_flag."
    exit 1
fi
if ! [[ -v pattern ]]; then
    >&2 echo "Missing required argument $pattern_flag."
    exit 1
fi
if ! [[ -v evalue ]]; then
    >&2 echo "Missing required argument $evalue_flag."
    exit 1
fi
if ! [[ -v sample_regex ]]; then
    >&2 echo "Missing required argument $sample_regex_flag."
    exit 1
fi
if ! [[ -v gene_regex ]]; then
    >&2 echo "Missing required argument $gene_regex_flag."
    exit 1
fi
if ! [[ -v keep_all ]]; then
    >&2 echo "Missing required argument $keep_all_flag."
    exit 1
fi
if ! [[ -v output_graph ]]; then
    >&2 echo "Missing required argument $output_graph_flag."
    exit 1
fi
cd_flag=
if [[ -v cache_dir ]]; then
    cd_flag="-D $cache_dir"
fi
if [ "$keep_all" = true ]; then
    ka_flag="-k"
else
    ka_flag=
fi
set -e
set -x
# Get top genes.
bash select_top_genes/select_top_sets_all.sh -t "$transcripts_fn" -n "$top_n" \
     -o "$out_dir_1" -p "$pattern" -j "$jobs" "${dirlist[@]}"
# Get gene matches tables.
python find_all_pairs.py -i "$out_dir_1"/*.fasta -O "$out_dir_2" \
       -r "$gene_regex" -R "$sample_regex" -e "$evalue" -n "$top_N" $ka_flag \
       $cd_flag -j "$jobs"
# Get gene matches graph.
shopt -s nullglob
python build_graph.py -i "$out_dir_2"/*.{pkl,h5} -o "$output_graph"
if [[ -v clean_level ]]; then
    if [ "$clean_level" -ge 1 ]; then
	rm -r "$out_dir_1"/*.fasta
	if [ -n "$cache_dir" ]; then
	    rm -r "$cache_dir"
	fi
    fi
    if [ "$clean_level" -ge 2 ]; then
	rm -r "$out_dir_2"/*.{pkl,h5}
    fi
fi
echo "$output_graph"
echo "${#dirlist[@]}"
