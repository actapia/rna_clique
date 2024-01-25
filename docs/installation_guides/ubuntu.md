# Installing RNA-clique for Ubuntu

This guide will walk you through installing RNA-clique on Ubuntu. The guide
assumes you have some basic familiarity working with Ubuntu using the command
line. The guide is designed for Ubuntu 23.10 but should also work on other
recent versions of Ubuntu, including the latest LTS release as of this writing,
Ubuntu 22.04.

## Download RNA-clique

For now, you can download RNA-clique by cloning its GitHub repository.

```bash
git clone --recurse-submodules https://github.com/actapia/rna_clique
```


## Installing BLAST

You can install a somewhat recent version of NCBI BLAST+ through APT.

```bash
sudo apt install ncbi-blast+ ncbi-blast+ ncbi-blast+-legacy
```

## Installing CPAN/cpanm

The BioPerl Bio::SeqIO module isn't present in Ubuntu default repositories, so
we will install `cpanm` to help us install that library.

```bash
sudo apt install cpanm
```

## Installing Bio::SeqIO

Now that `cpanm` is installed, we install Bio::SeqIO with the following command:

```bash
sudo cpanm install Bio::SeqIO
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

When asked about initializing Miniconda for the current shell, say yes.

When the installation has finished, you will need to restart your shell.

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
