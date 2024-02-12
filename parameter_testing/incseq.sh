#!/usr/bin/env bash
seq "$@"
if [ "$#" -gt 2 ] && [ $(( ("$3" - "$1") % "$2" )) -ne 0 ]; then
    echo "$3"
fi
