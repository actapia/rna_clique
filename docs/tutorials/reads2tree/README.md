# From RNA-seq reads to a phylogenetic tree with RNA-clique

This tutorial explains how to obtain a distance matrix and phylogenetic tree for
a set of samples for which we have RNA-seq reads. The purpose of this guide is
to show a typical workflow using RNA-clique.

The tutorial assumes that RNA-clique has already been downloaded and that its
dependencies and conda environment have already been installed.

## Setup

In addition to RNA-clique, this tutorial requires the following software:

* [sratoolkit](https://github.com/ncbi/sra-tools)
* [download_sra](https://github.com/actapia/download_sra)
* [SPAdes](https://github.com/ablab/spades)

This section provides brief installation instructions for each piece of
software. More detailed instructions may be found at each program's GitHub
repository.

This tutorial also expects the reader to have a POSIX-compatible shell like
`bash` or `zsh`, and common command-line utilities like `wget` or `curl`,
`tail`, `cut`, and `basename`. Some commands can also be run via GNU Parallel,
but this is optional.

It is recommended that the software be downloaded somewhere outside the
RNA-clique git repository. For example, you may wish to put the software in your
`~/Documents` directory.

### sratoolkit

Download the appropriate `sratoolkit` binaries for your system.

=== "Ubuntu"
    ```bash
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
    ```
=== "macOS (Intel)"
    ```zsh
    curl -L -O https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz
    ```
=== "macOS (Apple Silicon)"
    ```
    curl -L -O https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac-arm64.tar.gz
    ```

Then, extract the downloaded tar file.

```bash
tar xzvf sratoolkit.current-*.tar.gz
```

Add the `bin` directory of the extracted archive to your `PATH`.

```bash
export PATH="$PATH:$(realpath sratoolkit*/bin)"
```

### download_sra

!!! note
    If you downloaded this software from Zenodo, the test data SRA files are
    included in the software release zip. If you don't want to download them
    yourself, you don't need `download_sra` and can skip this step.

Make sure the `rna-clique` Conda environment created during the RNA-clique
installation is activated.

```bash
conda activate rna-clique
```

Install dependencies for `download_sra`.

```bash
conda install lxml requests
```

Clone the `download_sra` Git repository.

```bash
git clone https://github.com/actapia/download_sra
```

Add the repository root to your `PATH`.

```bash
export PATH="$PATH:$PWD/download_sra"
```

### SPAdes

!!! note
    If you downloaded this software from Zenodo, the assembled transcriptomes
    are included in the release zip. If you don't want to assemble the
    transcriptomes yourself, you can skip downloading SPAdes.

Download the SPAdes assembler. As of this writing, 4.0.0 is the newest
version.

=== "Ubuntu"
    ```bash
    wget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz
    ```
=== "macOS"
    ```zsh
    curl -L -O https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Darwin-$(uname -m).tar.gz
    ```

Extract the archive:

```bash
tar xzvf SPAdes-*.tar.gz
```

Add the SPAdes `bin` directory to your `PATH`.

```bash
export PATH="$PATH:$(realpath SPAdes-*/bin)"
```

### GNU Parallel (optional)

We can use GNU Parallel to run multiple SPAdes jobs
simultaneously. Parallelization can speed up these parts on systems with more
than one logical core ("thread").

=== "Ubuntu"
    ```bash
	sudo apt install parallel
	```
=== "macOS"
	```zsh
	brew install parallel
	```


## Creating a directory for our work

We will organize our files for this tutorial nicely into a single directory, but
first, we want to keep track of the location of RNA-clique. For this tutorial,
we will assume that RNA-clique is located at the path specified in the
`RNA_CLIQUE` environment variable. To set this variable easily, you can `cd`
into the Git repository and run

```bash
export RNA_CLIQUE=$PWD
```

We recommend putting the directory for this tutorial *outside* of the RNA-clique
Git repository. For the rest of this tutorial, we will assume that the tutorial
directory is at the path in the `TUTORIAL_DIR` environment variable. If you are
in the root of the Git repository, you could run 

```bash
cd ..
mkdir tutorial
cd tutorial
export TUTORIAL_DIR=$PWD
```


## Obtaining sequence data

!!! note
    If you downloaded this software from Zenodo, you already have the SRA files
    in the `test_data/sra` of the repository. Instead of completing this step,
    you can extract the provided data by running
	
	```bash
	for f in "$RNA_CLIQUE/test_data/sra/*"; do
	    fasterq-dump "$f"
	done
	```

!!! note
    Before proceeding, check that your environment variables are set to the 
    correct values!

    ```bash
    echo "$RNA_CLIQUE"
    echo "$TUTORIAL_DIR"
    ```

    If either of these does not print a path, you have forgotten to set the
    appropriate environment variables.

We will start out with RNA-seq reads for six samples of tall fescue (*Lolium
arundinaceum*). These samples are a subset of the sixteen used in the paper
*RNA-clique: A method for computing genetic distances from RNA-seq data*. The
sample metadata is shown below.

| SRA Accession                                                                                 | Genotype |
|-----------------------------------------------------------------------------------------------|----------|
| [SRR2321388](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browserSRR2321388acc=SRR2321388) | CTE46    |
| [SRR2321385](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browserSRR2321385acc=SRR2321385) | CTE46    |
| [SRR8003761](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browserSRR8003761acc=SRR8003761) | CTE27    |
| [SRR8003762](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browserSRR8003762acc=SRR8003762) | CTE27    |
| [SRR7990321](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browserSRR7990321acc=SRR7990321) | FATG4    |
| [SRR8003736](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browserSRR8003736acc=SRR8003736) | NTE      |

The list of accessions is in the file
[`tall_fescue_accs.csv`](./tall_fescue_accs.csv). We can download these all
easily using the `download_sra` tool.

First, change to your `TUTORIAL_DIR`.

```bash
cd "$TUTORIAL_DIR"
```

Then, run `download_sra.sh` on the tall fescue accessions. We will use the `-j`
flag to run the download with multiple jobs and the `-r` flag to remove the SRA
files after downloading and extracting.

```bash
tail -n+2 "$RNA_CLIQUE/docs/tutorials/reads2tree/tall_fescue_accs.csv" | \
    cut -d, -f1 | download_sra.sh -j 0 -r 
```

Verify that the FASTQ files have been extracted.

```bash
ls SRR*.fastq
```

## Assembling transcriptomes

!!! note
    If you downloaded this software from Zenodo, you already have the assembled
    transcriptomes in the `test_data/assemblies/` directory in the root of the
    repository and can simply move them instead of running SPAdes
	
	```bash
	mv "$RNA_CLIQUE/test_data/assemblies" out
	```
	
!!! note
    Assembling the transcriptomes requires at least 16 GB of memory. If you have
    insufficient memory or otherwise need to skip this step, you can download
    the assemblies instead:
	
	```bash
	wget "http://rna-clique-data.s3-website.us-east-2.amazonaws.com/transcripts.zip"
	mkdir out
	unzip transcripts.zip -d out
	```

Ordinarily, we would need a quality control step before proceeding to assembly,
but we will skip that for this tutorial.

We will use the rnaSPAdes mode of the SPAdes assembler to assemble the reads we
just downloaded into transcriptomes for each of the six samples. 

If your computer has sufficient resources, we will perform these assemblies in
parallel to save time. We estimate we need 16 GB of memory to assemble one
of these transcriptomes with three threads, so we want to run with no more than
$&lfloor; m / 16 &rfloor;$ jobs, where $m$ is the memory your computer has, in
GB. 

Determining exactly how many jobs and how many threads are optimal requires
trial and error; you may need to simply guess how many are appropriate for your
computer and retry if you run out of memory.

On a computer with over 120 GB of memory, we can run 6 jobs with 3 threads
safely.

=== "With parallel"
    ```bash
    parallel --jobs 6 spades.py --rna -o out/{/.} -s {} -t 3 -m 120 ::: *.fastq
    ```
=== "Without parallel"
    ```bash
    for f in *.fastq; do
        b="$(basename "$f")"; 
        fn="${b%%.*}";
        spades.py --rna -o "out/$fn" -s "$f" -t 3 -m 120;
    done
    ```

The assembled transcriptomes will be located at `transcripts.fasta` in
directories corresponding to their samples names under the `out` directory.

## Running phase 1 of RNA-clique

First, return to the RNA-clique repository root.

```bash
cd "$RNA_CLIQUE"
```

Then, run `filtering_step.py` on the transcriptomes we just assembled. This
script will perform "phase 1" of RNA-clique, which involves:

1. Selecting the top $n$ genes for each transcriptome
2. BLASTing each sample's top genes against every other's to get gene matches
   tables
3. Building the gene matches graph from the gene matches tables

<!-- We need to create a directory in which RNA-clique can write its output. We'll -->
<!-- create a directory named `rna_clique_out` in the `$TUTORIAL_DIR`. -->

<!-- ```bash -->
<!-- mkdir "$TUTORIAL_DIR/rna_clique_out" -->
<!-- ``` -->

Previous tests with this data revealed that $n = 50000$ is a good setting, so we
will use that value.

```bash
python filtering_step.py -O "$TUTORIAL_DIR"/rna_clique_out \
                         -n 50000 \
                         "$TUTORIAL_DIR"/out/*
```

Verify that the `graph.pkl` file was created in the output directory.

```bash
ls "$TUTORIAL_DIR/rna_clique_out/graph.pkl"
```

## Results

### Getting a tree

If you want a tree, you can create one using RNA-clique and Biopython. The code
below, also found in `docs/tutorials/reads2tree/make_tree.py`, computes the
distance matrix from the `graph.pkl` and `od2/*.h5` (or `od2/*.pkl`) files and
constructs a tree using the neighbor-joining algorithm. The tree is also rooted
at its midpoint. The tree is saved to `nj_tree.tree`, and a visualization is
saved to `nj_tree.svg` in the `rna_clique_out` directory.

```python
--8<-- "docs/tutorials/reads2tree/make_tree.py"
```


The script requires some modules found in the root of the RNA-clique repository,
so you can run it as follows:

```bash
PYTHONPATH='.' python docs/tutorials/reads2tree/make_tree.py
```

### Getting a PCoA plot

We can use `scikit-bio` to create a PCoA plot from our distance matrix. To
distinguish points by genotype, we will need to use the metadata for the samples
stored at `$RNA_CLIQUE/docs/tutorials/reads2tree/tall_fescue_accs.csv`.

The code below draws a 3D and 2D PCoA plot and stores the results as SVG files
in the `rna_clique_out` directory as `pcoa_3d.svg` and `pcoa_2d.svg`,
respectively. The code can also be found at
`docs/tutorials/reads2tree/make_pcoa.py`.

```python
--8<-- "docs/tutorials/reads2tree/make_pcoa.py"
```

The example can be run as follows from the root of the RNA-clique repository.

```bash
PYTHONPATH="." python docs/tutorials/reads2tree/make_pcoa.py
```

### Getting a heatmap

We can use the `draw_heatmap` function of RNA-clique to display a similarity or
distance matrix as a heatmap. The function uses the Seaborn `heatmap` function
behind the scenes, and arbitrary arguments given to `draw_heatmap` will be
passed to Seaborn.

The code below is also found in `docs/tutorials/reads2tree/make_heatmap.py`. It
draws a heatmap and saves the resulting figure in the `rna_clique_out` directory
as `distance_heatmap.svg`. 

```py
--8<-- "docs/tutorials/reads2tree/make_heatmap.py"
```


To generate a heatmap using this code, you can run the Python script as follows
from the RNA-clique repository root.

```bash
PYTHONPATH="." python docs/tutorials/reads2tree/make_heatmap.py
```
