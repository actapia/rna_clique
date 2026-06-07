#!/usr/bin/env zsh
branch="main"
if [ "$#" -gt 2 ]; then
    >&2 echo "This script accepts at most one argument."
    exit 1
fi
if [ "$#" -eq 1 ]; then
    branch="$1"
fi
set -e
cd
CI=1 /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
eval "$(/opt/homebrew/bin/brew shellenv zsh)"
brew install blast python@3.14
python3.14 -m venv rna_clique_venv
. rna_clique_venv/bin/activate
git clone -b "$branch" --recurse-submodules https://github.com/actapia/rna_clique
cd rna_clique
python -m pip install build
python -m build
python -m install dist/*.whl
bash tests/verify_install/test_install.sh
