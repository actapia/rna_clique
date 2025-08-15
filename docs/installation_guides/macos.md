# Installing RNA-clique for macOS

This guide will walk you through installing RNA-clique on macOS. Since
RNA-clique currently has no graphical interface, basic command-line skills are
needed for both installation and use of the software.

This guide is designed for macOS Sequoia 15 on Intel Mac, but the steps in this
guide have also been tested on macOS Ventura 13 on Intel and macOS Sonoma 14 on
Apple Silicon via GitHub Actions runners. The guide likely works on other recent
versions of macOS, but such configurations are untested.

## Installing command line developer tools

Installation requires some commands (e.g., `git`) that are part of the macOS 
command line developer tools, which are not installed by default. To install
them, run the following in the Terminal.

```zsh
xcode-select --install
```

If the command line developer tools have not been installed already, you should
see a graphical dialog prompting you to confirm installation. Click "Install" 
and agree to the license agreement to begin installation .Wait for the install
to complete before continuing.

## Download RNA-clique

You can download RNA-clique by cloning its GitHub repository. (Alternatively,
you can extract the release zip or tarball, but be mindful that the repository
root will have a different name than the one presented here.)

<!--{{clone_command(git_branch()) | code_fence("zsh") | comment_surround}}{{empty("-->
```bash
git clone --recurse-submodules https://github.com/actapia/rna_clique
```
<!--")}}-->

## Installing Homebrew

Homebrew (`brew`) is a package manager for macOS---it allows the user to install
and update software packages via the command line. We will use `brew` to install
some of the dependencies of RNA-clique, but we need to install the package
manager itself first.

```zsh
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

You will most likely need to enter your password to install Homebrew as root.
The install script will also ask you to confirm that you want to install the
needed software.

## Installing Bash

We will install a few packages using Homebrew.

First, we will install **Bash**, which is a Unix shell. The version of Bash that
comes installed on macOS is very out of date, but RNA-clique requires some new
features of Bash.

Second, we will install **NCBI BLAST+**, a popular implementation of the BLAST
local sequence alignment software.

Third, we will install **cpanminus**`, a Perl package manager that will help us
install the Perl dependencies.

```zsh
brew install bash blast cpanminus
```

## Installing GNU Parallel (optional)

Some parts of RNA-clique can use GNU Parallel to run multiple jobs 
simultaneously. Parallelization can speed up these parts on systems with more
than one logical core ("thread"). GNU parallel be installed via `brew`.

```zsh
brew install parallel
```

## Installing Bio::SeqIO

We can use `cpan` to install the Bio::SeqIO Perl module.

```zsh
sudo cpanm Bio::SeqIO
```

macOS may ask you to permit Terminal to access certain files and settings
on your computer. You can safely allow Terminal that access.

## Installing Array::Heap

Again, we can use `cpan`.

```zsh
sudo cpanm Array::Heap
```

## Installing Miniconda

The conda package and environment manager will be useful for setting up the
Python dependencies of RNA-clique. If you already have a conda installation, you
can skip to the next step, [Installing the rna-clique conda
environment](#installing-the-rna-clique-conda-environment) First, download the
Miniconda installer:

```zsh
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-$(uname -m).sh
```

Then, run the installer.

```zsh
bash Miniconda3-latest-MacOSX-$(uname -m).sh
```

Agree to the license terms. Unless you have a reason to install somewhere else,
you can install to the default location Miniconda suggests.

When asked about initializing Miniconda for the current shell, say yes.

When the installation has finished, you will need to restart your shell.

### Installing the rna-clique conda environment

Make sure you are in the `rna_clique` repository that you cloned with `git`.

```zsh
cd rna_clique
```

Then, create a `conda` environment from the `environment.yml` file in the root
of the repository.

```zsh
conda env create -f environment.yml --name rna-clique
```

When asked to confirm that you want to install the dependencies, say yes.

When creation of the new environment has finished, you can activate the
environment to ensure you are running Python with the needed dependencies:

```zsh
conda activate rna-clique
```
