# RNA-clique

[![DOI](https://zenodo.org/badge/747462438.svg)](https://zenodo.org/badge/latestdoi/747462438)

This is the repository for RNA-clique, a tool for computing pairwise genetic
distances from RNA-seq data. The software accepts as input assembled
transcriptomes from two or more samples and produces as its output a matrix
containing pairwise distances ranging from 0 to 1.

## Installation

This software is written in Python, Perl, and Bash. The software additionally
requires NCBI BLAST+ and several Python and Perl libraries. The section below
lists the software requirements, and [guides](#installation-guides) are provided
for installation of specific systems.

### Requirements

The requirements below represent one tested configuration. This software may
work with different versions of the dependencies listed, but such configurations
are considered untested.

#### Main software

* Python 3.12
* Perl 5.38.0
* Bash 5.2.32
* ncbi-blast 2.12.0+
* Python libraries
  * tqdm 4.66.5
  * pandas 2.2.2
  * pytables 3.10.1
  * joblib 1.4.2
  * networkx 3.3
  * [simple-blast](https://github.com/actapia/simple_blast) 0.0.5
* Perl libraries
  * Bio::SeqIO 1.7.8
    
#### Phylogenetics and Visualization

* Python libraries
  * BioPython 1.84
  * matplotlib 3.9.2
  * seaborn 0.13.2
  * scikit-bio 0.6.2
  
#### Sequence simulation (needed for testing installation)

* Python libraries
  * more-itertools 10.4.0
  * PyYAML 6.0.2
  * SciPy 1.14.1
  

### Installation guides

* [Ubuntu](https://actapia.github.io/rna_clique/installation_guides/ubuntu)
* [macOS](https://actapia.github.io/rna_clique/installation_guides/macos)

### Testing the installation

After you have downloaded RNA-clique and installed its dependencies, you can
test your installation using a script included in this repository.

First, make sure your rna-clique conda environment is active.

```bash
conda activate rna-clique
```

Then, if you are in the root the repository, you can run the following command
to begin the test script.

```bash
bash tests/verify_install/test_install.sh && echo "Success!"
```

The script generates a small test dataset and runs RNA-clique on the generated
data. On a modern desktop with one thread, the test should take around one
minute to complete. On machines with multiple threads, the test script should
take advantage of parallelism to complete the test more quickly.

If you ran the script with the above command and see "Success!," then the
installation was succesful. Otherwise, you will need to investigate the output
of the test script to see what failed and why.

If the test script fails despite having a correct installation, you should
submit a bug report on GitHub at https://github.com/actapia/rna_clique/issues .

## Command-line usage

Running RNA-clique broadly involves two phases. In the first phase, the
transcriptomes are aligned, and the gene matches graph is built. In the second
phase, the gene matches graph is used to filter the BLAST alignments, and the
pairwise distances are calculated using the filtered alignment statistics.

Although there are only two phases, each phase may be performed by multiple
scripts. To simplify usage of this program, we have provided a script,
`typical_filtering_step.sh` that may be used to easily perform the first phase
with typical parameter settings.

### Phase 1: Building the gene matches graph

The simplest usage of `typical_filtering_step.sh` provides only an output
directory, a value for the number of top genes to select (`n`), and the
directories containing the transcriptomes to be analyzed.

The script assumes that the transcriptomes are stored in FASTA files with the
*identical names* in different directories. By default, RNA-clique assumes the
files are all named `transcripts.fasta`, since this is the default output name
for the SPAdes assembler, but this behavior may be overridden by a command-line
argument. 

For example, the first transcriptome might be located at
`sample1/transcripts.fasta`, and the second might be located
at`sample2/assembly.fasta`.

**WARNING: The transcriptomes are identified internally by the names of the
directories in which they are contained, so every transcriptome must be located
in a directory with a unique name!**

```bash
bash typical_filtering_step.sh -n TOP_GENES -o OUTPUT_DIR DIR1 DIR2 ...
```

In the command above, `TOP_GENES` must be replaced by the number of top genes to
select. `OUTPUT_DIR` should be replaced by the path to the directory in which to
store output files. `DIR1 DIR2 ...` are the paths to the directories containing
the transcriptomes to be analyzed. The transcriptomes are assumed be located at
files named `transcripts.fasta` within those directories.

When the script finishes, it creates `graph.pkl` in the specified output
directory. `graph.pkl` is a Python pickle file representing the constructed
gene matches graph.

The script also stores HDF5 files (formerly Python pickles) for the pairwise
BLAST results. The BLAST results can be found in the `od2` subdirectory of the
output directory.

### Phase 2: Calculating distances

The `filtered_distance.py` Python script may be used to compute distances or
similarities from a gene matches graph. Basic usage of the command requires
only that we provide the pickles for the gene matches graph and the HDF5 files
for the pairwise BLAST results.

```bash
python filtered_distance.py -g GRAPH -c COMPARISONS_DIR/*.h5
```

In the above command, GRAPH should be the path to the `graph.pkl` created in the
first phase, and COMPARISONS_DIR should be the directory that contains the BLAST
result HDF5 files. (This will be the `od2` subdirectory of the output directory
from Phase 1 if you used the `typical_filtering_step.sh` script.)

The script outputs a genetic similarity matrix to standard output by default. To
get a distance matrix, you can provide the `-o dis` option to
`filtered_distance.py`.

Both the rows and the columns of the matrix are sorted alphabetically by sample
ID. To print the order of the samples to standard error, you can provide the
`-l` flag to `filtered_distance.py`.

### Downstream analyses

The `filtered_distance.py` script prints the calculated matrix to the standard
output, so you can use redirection or pipes to save the results to a file. You
could then use the matrix in any downstream application capable of loading
arbitrary matrices from files.

For example, if you output the matrix to a file named `distances`, you could
load the matrix in R using the following code:

```R
dis <- as.matrix(read.table("distances", sep=" "))
```

If you are using Python, you may wish to skip Phase 2 and [use RNA-clique
directly in your code instead](#using-rna-clique-in-python-code).


## Using RNA-clique in Python code

Since parts of Phase 1 are implemented in Bash and Perl, there is currently no
official way to perform Phase 1 from custom Python code, but since Phase 2 is
written exclusively in Python, we describe an official way of performing that
phase in custom code here.

The `filtered_distance.py` is a straightforward command-line interface to
RNA-clique's `SampleSimilarity` class. If you are planning to use the distance
matrix computed by RNA-clique in a downstream analysis implemented in Python, it
may be easier to simply use `SampleSimilarity` directly instead of running
`filtered_distance.py` and loading the calculated matrix back into Python.

The `SampleSimilarity` class is ordinarily constructed with the following
arguments:
 
 * A NetworkX graph representing the gene matches graph.
 * An Iterable of pairs to be interpreted as a mapping between unordered pairs
   of sample names and dataframes containing the parsed BLAST results for the
   pairs' comparisons. (For example, if you had a `dict` `d` mapping
   `frozenset`s containing pairs of sample names to the BLAST results for those
   pairs of samples, then `d.items()` would be an appropriate actual parameter.)
 * An optional integer indicating how many samples are present. (If no value is
   provided for this parameter, it will be calculated automatically.)
   
This constructor may be difficult to use if you are simply given the file paths
to the pickles for the graph and the BLAST comparisons, so `SampleSimilarity`
also includes a `from_filenames` classmethod that constructs a
`SampleSimilarity` object using the following parameters:

* A path to the pickle containing the gene matches graph.
* A list of paths to the BLAST comparisons.
* An optional `bool` indicating whether to keep the comparison dataframes after
  the distance computation has finished.

`SampleSimilarity` computes the similarities lazily; it won't do any of the
time-consuming work until it's asked for results. `SampleSimilarity` offers
several ways to view the results.

`get_similarities` and `get_dissimilarities` return a `MultisetKeyDict` that
maps unordered pairs of sample names to their similarities or dissimilarities
(distances), respectively. The `MultisetKeyDict` for the  similarities may also
be accessed via the `similarities` property, but note that there is no
corresponding `dissimilarities` or `distances` property.

`get_similarity_matrix` and `get_dissimilarity_matrix` offer the same
information as matrices&mdash;i.e., NumPy arrays. The returned matrices are
symmetric, but the distances are only computed in one direction by
RNA-clique. The rows and columns of the matrix correspond to the samples sorted
alphabetically by name. To get this ordered list of samples, you can use the
`samples` property of `SampleSimilarity`.

## License

All code except tests/test_install/macos.sh is licensed under the MIT license,
which may be found at LICENSE.MIT at the root of this repository.

[tests/test_install/macos.sh](https://github.com/actapia/rna_clique/blob/main/tests/test_install/macos.sh)
is licensed under the Creative Commons Attribution-ShareAlike 4.0 License, which
may be found at LICENSE.CC-BY-SA-4.0 at the root of this repository.

A machine-readable copyright file in Debian format may also be found at
[copyright](https://github.com/actapia/rna_clique/blob/main/copyright).

## Additional documentation

* [Command-line usage guide](https://actapia.github.io/rna_clique/usage)
* [Tutorial: From RNA-seq reads to a phylogenetic tree with RNA-clique](https://actapia.github.io/rna_clique/tutorials/reads2tree)
