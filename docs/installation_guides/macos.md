# Installing RNA-clique for macOS

This guide will walk you through installing RNA-clique on macOS. Since
RNA-clique currently has no graphical interface, basic command-line skills are
needed for both installation and use of the software.

This guide is designed for macOS 26 Tahoe and has also been tested on macOS 15
Sequoia via GitHub Actions runners. Newer versions of macOS will likely work,
but older versions probably will not.

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

After the install has finished, run the commands that the installer suggests to
add Homebrew to your `PATH`, and the restart your shell.

## Installing BLAST and Python

We will use `brew` to install **NCBI BLAST+**, a popular implementation of the
BLAST local sequence alignment software.

We will also install a new version of Python since macOS comes pre-installed
with a version that is too old to work with RNA-clique.

```zsh
brew install blast python@3.13
```

## Creating a virutal environment

It is good practice to install our software to a separate Python virtual
environment. Software installed in the virtual environment will only be usable
when we have the virtual environment "activated," so it won't clutter or
conflict with the other software installed on the system. 

```zsh
python3.13 -m venv rna_clique_env
```

Then, activate the environment with the virtual environment's `activate` script.

```zsh
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
