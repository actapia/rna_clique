# Building the RNA-clique Python package

As of this writing, RNA-clique is not yet available from the Python Package
Index repositories, so it must be built and installed manually. This document
describes the steps necessary to build an RNA-clique wheel and install it using
`pip`.

## Ensure you are in the correct environment

First, make sure you are in the environment in which you want to install
RNA-clique. If you followed the installation guide, you should have a Python
virtual environment you can use.

## Downloading the RNA-clique source code

If you haven't already, download RNA-clique. There are a few options for
downloading RNA-clique, including downloading archives from Zenodo (full
releases only) or GitHub. This document, however, will explain how to download
RNA-clique via `git`. If you use a different method, note that the name of the
directory containing RNA-clique might differ from the name assumed in this
tutorial. 

!!! tip 
    If you've already installed `brew` on macOS, you should have `git`
    already. If you are on Ubuntu, you will need to install `git` via `apt`.
	
	```bash
	sudo apt update && sudo apt install git
	```
	
Choose a directory in which to download the `rna_clique` source code, then run

<!--{{clone_command(git_branch()) | code_fence("zsh") | comment_surround}}{{empty("-->
```bash
git clone --recurse-submodules https://github.com/actapia/rna_clique
```
<!--")}}-->

## Install build

We will use the `build` module to build RNA-clique, but we need to install that
module first.

```bash
python -m pip install build
```

## Build RNA-clique

First, change to the directory containing the RNA-clique source code.

```bash
cd rna_clique
```

Then, run

```bash
python -m build
```

## Installing the built wheel

Once `build` finishes, you can install the built wheel with

```bash
python -m pip install dist/*.whl
```
