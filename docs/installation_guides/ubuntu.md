# Installing RNA-clique for Ubuntu

This guide will walk you through installing RNA-clique on Ubuntu. The guide
assumes you have some basic familiarity working with Ubuntu using the command
line. The guide is designed for Ubuntu 23.10 and has also been tested on Ubuntu
22.04 and 24.04. The guide likely works on other recent versions of Ubuntu as
well, but such configurations are untested.

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


```bash
sudo apt install git wget ncbi-blast+ ncbi-blast+ ncbi-blast+-legacy
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

You can download RNA-clique by cloning its GitHub repository. (Alternatively,
you can extract the release zip or tarball, but be mindful that the repository
root will have a different name than the one presented here.)

<!--{{clone_command(git_branch()) | code_fence("bash") | comment_surround}}{{empty("-->
```bash
git clone --recurse-submodules https://github.com/actapia/rna_clique
```
<!--")}}-->

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
