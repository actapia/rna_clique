#!/usr/bin/env zsh
set -x
branch="main"
if [ "$#" -gt 2 ]; then
    >&2 echo "This script accepts at most one argument."
    exit 1
fi
if [ "$#" -eq 1 ]; then
    branch="$1"
fi

# Install Command Line Tools.
# Based on an answer from keen on Ask Different.
# https://apple.stackexchange.com/a/195963
readonly RETRIES=5
touch /tmp/.com.apple.dt.CommandLineTools.installondemand.in-progress
res=1
i=0
PROD="$(softwareupdate -l |
  grep "\*.*Command Line" |
  head -n 1 | awk -F"*" '{print $2}' |
  sed -e 's/^ *//' |
  tr -d '\n' | sed -e 's/^Label: //')"
if [ $? -ne 0 ]; then
    >&2 echo "Could not get Command Line Tools version from update server."
fi
while [ "$res" -ne 0 ] && [ "$i" -lt "$RETRIES" ]; do
    softwareupdate -i "$PROD" --verbose
    res=$?
    i=$((i+1))
done
if [ "$res" -ne 0 ]; then
    >&2 echo "Could not install command line tools. (Maybe retry?)"
    exit 1
fi
rm /tmp/.com.apple.dt.CommandLineTools.installondemand.in-progress

set -e
git clone -b "$branch" --recurse-submodules https://github.com/actapia/rna_clique
if [ "$PASSWORDLESS" -ne 1 ] && which sudo && ! sudo -n -v; then
    echo "$password" | sudo -S -v
fi
CI=1 /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
brew install bash blast
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-$(uname -m).sh
bash Miniconda3-latest-MacOSX-$(uname -m).sh -b -f
"$HOME/miniconda3/bin/conda" init zsh
"$HOME/miniconda3/bin/conda" init bash
cd rna_clique
eval "$("$HOME/miniconda3/bin/conda" shell.zsh hook)"
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
conda env create -y -f environment.yml --name rna-clique
conda activate rna-clique
bash tests/verify_install/test_install.sh
