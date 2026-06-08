# Installing RNA-clique for Ubuntu

This guide will walk you through installing RNA-clique on Ubuntu. The guide
assumes you have some basic familiarity working with Ubuntu using the command
line. The guide is designed for Ubuntu 26.04 and has also been tested on Ubuntu
24.04 via GitHub Actions runners.

## Installing dependencies via APT

RNA-clique requires that NCBI BLAST+, a popular implementation of the BLAST
sequence alignment algorithm, be installed. 

For Python, we also need to be able to create virtual environments, so we will
need to install `venv`. Some Python packages we will install also need a C++
compiler, so we will install the common `g++` compiler. We can install both
using Ubuntu's default package manager, APT.

First, update the package lists.

```bash
sudo apt update
```

Then, install NCBI BLAST+ and `g++`.

```bash
sudo apt install ncbi-blast+ ncbi-blast+-legacy g++ python3-venv
```

## Creating a virutal environment

It is good practice to install our software to a separate Python virtual
environment. Software installed in the virtual environment will only be usable
when we have the virtual environment "activated," so it won't clutter or
conflict with the other software installed on the system. 

```bash
python3 -m venv rna_clique_env
```

Then, activate the environment with the virtual environment's `activate` script.

```bash
. rna_clique_env/bin/activate
```

You should see that `(rna_clique_env)` appears at the beginning of your prompt,
letting you know.

## Installing RNA-clique

You can install RNA-clique from the Python Package Index using
`pip`. Alternatively, you can [build the package yourself](../building.md) from
the Git repository.

```python
python -m pip install rna-clique
```

!!! note
    When you log out and log back into your shell, you will need to reactivate
    the virtual environment using the `activate` script.
