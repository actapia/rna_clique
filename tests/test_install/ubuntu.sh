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
sudo NEEDRESTART_MODE=a apt install -y git wget ncbi-blast+ ncbi-blast+ ncbi-blast+-legacy \
     build-essential cpanminus
sudo NEEDRESTART_MODE=a apt install -y bioperl --no-install-recommends
sudo cpanm Array::Heap
wget --no-verbose https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh -b -f
"$HOME/miniconda3/bin/conda" init bash
eval "$("$HOME/miniconda3/bin/conda" shell.bash hook)"
git clone -b "$branch" --recurse-submodules https://github.com/actapia/rna_clique
cd rna_clique
conda env create -y -f environment.yml --name rna-clique
conda activate rna-clique
bash tests/verify_install/test_install.sh
