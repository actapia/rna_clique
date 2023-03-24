#!/usr/bin/env bash
readarray -t params
i=0
while IFS= read -u 5 -r line; do
    echo "$line ${params[$i]}"
    i=$((i+1))
done 5< <(printf '%s\n' "${params[@]}" | python gen_params.py "$@" | python add_output_params.py -i | parallel --colsep " " "-j$TRY_PARAMS_JOBS" bash stats_for_params.sh)
