# RNA-clique

This is the repository for RNA-clique, a tool for computing pairwise genetic distances from RNA-seq data. The software accepts as input assembled transcriptomes from two or more samples and produces as its output a matrix containing pairwise distances ranging from 0 to 1.

## Installation

This software is written in Python, Perl, and Bash. The software additionally requires NCBI BLAST+ and several Python and Perl libraries. The section below lists the software requirements, and [guides](#installation-guides) are provided for installation of specific systems.

### Requirements

The requirements below represent one tested configuration. This software may work with different versions of the dependencies listed, but such configurations are considered untested.

#### Main software

* Python 3.11
* Perl 5.36.0
* Bash 5.2.15
* ncbi-blast 2.12.0+
* Python libraries
  * tqdm
  * pandas
  * joblib
  * networkx
  * [simple-blast](https://github.com/actapia/simple_blast)
* Perl libraries
  * Bio::SeqIO
  
  
#### Phylogenetics and Visualization

* Python libraries
  * BioPython
  * matplotlib
  * seaborn
  

### Installation guides

* [Ubuntu](docs/installation_guides/ubuntu.md)

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
bash tests/test_install/test_install.sh && echo "Success!"
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

## Usage

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

**WARNING: The transcriptomes are identified internally by their filenames, so
every transcriptome must have a unique filename!**

```bash
bash typical_filtering_step.sh -n TOP_GENES -o OUTPUT_DIR TRANSCRIPTOME1 TRANSCRIPTOME2 ...
```

In the command above, `TOP_GENES` must be replaced by the number of top genes to
select. `OUTPUT_DIR` should be replaced by the path to the directory in which to
store output files. `TRANSCRIPTOME1 TRANSCRIPTOME2 ...` are the paths to the
transcriptomes to be analyzed.

When the script finishes, it creates `graph.pkl` in the specified output
directory. `graph.pkl` is a Python pickle file representing the constructed
gene matches graph.

The script also creates Python pickles for the pairwise BLAST results. The BLAST
results can be found in the `od2` subdirectory of the output directory.

### Phase 2: Calculating distances

The `filtered_distance.py` Python script may be used to compute distances or
similarities from a gene matches graph. Basic usage of the command requires
only that we provide the pickles for the gene matches graph and the pairwise
BLAST results.

```bash
python filtered_distance.py -g GRAPH -c COMPARISONS_DIR/*.pkl
```

In the above command, GRAPH should be the path to the `graph.pkl` created in the
first phase, and COMPARISONS_DIR should be the directory that contains the BLAST
result pickles. (This will be the `od2` subdirectory of the output directory
from Phase 1 if you used the `typical_filtering_step.sh` script.)

The script outputs a genetic similarity matrix to standard output by default. To
get a distance matrix, you can provide the `-o dis` option to
`filtered_distance.py`.
