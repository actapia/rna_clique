#!/usr/bin/env bash
cat >&2 <<EOF
typical_filtering_step.sh is deprecated and is now a wrapper for
filtering_step.py. Please use filtering_step.py instead.
EOF
args=( "$@" )
for ((i=0;i<"${#args[@]}";i++)); do
    case "${args[i]}" in
	"-o")
	    echo >&2 "Interpreting -o as --output-dir (now -O)."
	    args[$i]="-O"
	    ;;
	*)
	    :
	    ;;
    esac
done
python filtering_step.py "${args[@]}"
