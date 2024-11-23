# Installing RNA-clique for Ubuntu

This guide will walk you through installing RNA-clique on Ubuntu. The guide
assumes you have some basic familiarity working with Ubuntu using the command
line. The guide is designed for Ubuntu 23.10 but should also work on other
recent versions of Ubuntu, including the latest LTS release as of this writing,
Ubuntu 22.04.

## Installing dependencies via APT

We will install some dependencies using Ubuntu's default package manager, APT.

First, update the package lists.

```bash
sudo apt update
```

Then, install the first set of packages. The table below describes what we will
install.

| Software          | Description                                                                            |
| ----------------- | -------------------------------------------------------------------------------------- |
| Git               | Version control system used by RNA-clique (usually installed by default)               |
| wget              | Command-line utility for downloading files from the web (usually installed by default) |
| NCBI BLAST+       | Popular implementation of the BLAST local sequence alignment algorithm                 |
| `build-essential` | Metapackage containing tools useful for building software , including Perl modules     |
| cpanminus         | Package manager for Perl                                                               |


```bash
sudo apt install git wget ncbi-blast+ ncbi-blast+ ncbi-blast+-legacy \
                 build-essential cpanminus
```

BioPerl can also be installed via APT. The package and its dependencies include
many recommendations that we don't need&mdash;you can specify
`--no-install-recommends` to avoid installing those.

```bash
sudo apt install bioperl --no-install-recommends
```

## Installing GNU Parallel (optional)

Some parts of RNA-clique can use GNU Parallel to run multiple jobs
simultaneously. Parallelization can speed up these parts on systems with more
than one logical core ("thread"). GNU parallel be installed via APT.

```bash
sudo apt install parallel
```

## Installing Array::Heap

We will use cpanmius to install the `Array::Heap` Perl  module, since it isn't
available in Ubuntu's default software repositories.

```bash
sudo cpanm Array::Heap
```

## Installing Miniconda

The conda package and environment manager will be useful for setting up the
Python dependencies of RNA-clique. If you already have a conda installation, you
can skip to the next step, [Installing the rna-clique conda
environment](#installing-the-rna-clique-conda-environment) First, download the
Miniconda installer:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Then, run the installer.

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

Unless you have a reason to install somewhere else, you can install to the
default locaton Miniconda suggests.

When asked about initializing Miniconda for the current shell, say yes.

When the installation has finished, you will need to restart your shell.

## Downloading RNA-clique

For now, you can download RNA-clique by cloning its GitHub repository.

```bash
git clone --recurse-submodules https://github.com/actapia/rna_clique
```

## Installing the rna-clique conda environment

Make sure you are in the `rna_clique` repository that you cloned with `git`.

```bash
cd rna_clique
```

Then, create a `conda` environment from the `environment.yml` file in the root
of the repository.

```bash
conda env create -f environment.yml --name rna-clique
```

When asked to confirm that you want to install the dependencies, say yes.

When creation of the new environment has finished, you can activate the
environment to ensure you are running Python with the needed dependencies:

```bash
conda activate rna-clique
```
