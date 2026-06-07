#!/usr/bin/env bash
branch="main"
if [ "$#" -gt 2 ]; then
    >&2 echo "This script accepts at most one argument."
    exit 1
fi
if [ "$#" -eq 1 ]; then
    branch="$1"
fi
if ! which sudo; then
    function sudo {
	"$@"
    }
fi
set -e
cd
sudo apt update
sudo NEEDRESTART_MODE=a apt install -y git ncbi-blast+ ncbi-blast+-legacy \
     g++
python -m venv rna_clique_venv
. venv/bin/activate
git clone -b "$branch" --recurse-submodules https://github.com/actapia/rna_clique
cd rna_clique
python -m pip install build
python -m build
python -m install dist/*.whl
bash tests/verify_install/test_install.sh
