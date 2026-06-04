# Requirements

This page aims to document all immediate software dependencies of
RNA-clique. Dependencies are divided by type in the following sections.

## Python dependencies

The 
{{file_link("`environment.yml`", "environment.yml")}} 
file provides a list of dependencies
written in a format usable by the conda package manager, but the format uses
YAML and is human-readable. The list is reproduced below.

```py
--8<-- "environment.yml"
```

## Other dependencies

RNA-clique requires [NCBI
BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
be installed. Many systems offer these executables through their package
managers.

[Testing the installation](testing_installation.md) requires Bash to run the
`test_install.sh` script.

Currently, a C++ compiler must also be installed for some Python extensions to
be built by `pip`. Building the extension modules on which RNA-clique depends
has been tested with g++ 15.2.0.
