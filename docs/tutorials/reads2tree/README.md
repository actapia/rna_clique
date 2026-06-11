# From RNA-seq reads to a phylogenetic tree with RNA-clique

This tutorial explains how to obtain a distance matrix and phylogenetic tree for
a set of samples for which we have RNA-seq reads. The purpose of this guide is
to show a typical workflow using RNA-clique.

The tutorial assumes that RNA-clique has already been installed.

The most resource-consuming part of this tutorial is downloading and assembling
the RNA-seq reads. If you would rather skip this step, you can instead download
pre-assembled transcriptomes by skipping to the
["Assembling transcriptomes"](#assembling-transcriptomes) section and following
the note under the section heading. If you downloaded this software from Zenodo,
these pre-assembled transcriptomes will already be included with the software.

## Background

The RNA-seq data used in this tutorial are from plants derived from the
"Kentucky 31" tall fescue (*Lolium arundinaceum*) cultivar and were originally
used in gene expression studies described in "[Transcriptome analysis and
differential expression in tall fescue harboring different endophyte strains in
response to water
deficit](https://doi.org/10.3835/plantgenome2018.09.0071)". Each set of RNA-seq
reads comes from a different individual, and although we will be using RNA-seq
reads for six individuals, the individuals have only four distinct
genotypes. Individuals with the same genotype are clones, and should thus have
almost no differences in their genomes.

Some of the individuals possess an endosymbiotic fungus, *Epichloë coenophiala*,
while others were treated to remove the fungus. The endophyte statuses of the
individuals, which are relevant for later tutorials such as the "[Quickly
computing subsets of existing analyses](../subsets/README.md)" tutorial, are
also shown in the metadata table below.

| SRA Accession                                                                                 | Genotype | Endophyte |
|-----------------------------------------------------------------------------------------------|----------|-----------|
| [SRR2321388](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browserSRR2321388acc=SRR2321388) | CTE46    | infected  |
| [SRR2321385](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browserSRR2321385acc=SRR2321385) | CTE46    | minus     |
| [SRR8003761](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browserSRR8003761acc=SRR8003761) | CTE27    | infected  |
| [SRR8003762](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browserSRR8003762acc=SRR8003762) | CTE27    | minus     |
| [SRR7990321](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browserSRR7990321acc=SRR7990321) | FATG4    | infected  |
| [SRR8003736](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browserSRR8003736acc=SRR8003736) | NTE      | infected  |

For the purpose of this tutorial, each individual (associated with a single set
of RNA-seq reads, and, eventually, a single assembled transcriptome) will be a
"sample," and the goal of the tutorial is to obtain a genetic distance matrix
that quantifies pairwise differences in the genomes of the samples. 

## Setup

This tutorial expects the reader to have a POSIX-compatible shell like `bash` or
`zsh`, and common command-line utilities like `tail`, `cut`, and
`basename`. Some commands can also be run via GNU Parallel, but this is
optional.

### Activate the environment

If you haven't already, activate the environment we created during the
setup.

```bash
. rna_clique_venv/bin/activate
```

If the environment has been activated, you should see `(rna_clique_venv)` appear
at the beginning of your prompt.

### Installing additional tools

!!! note
    If you are using pre-assembled transcriptomes, you can skip installing the
    software in this section.

In addition to RNA-clique, this tutorial requires the following software for
obtaining and assembling the sequence reads:

* [sratoolkit](https://github.com/ncbi/sra-tools)
* [download_sra](https://github.com/actapia/download_sra)
* [SPAdes](https://github.com/ablab/spades)

We also need `git` to download various repositories and `wget` or
`curl` to download software and data.

This section provides brief installation instructions for each piece of
software.

#### git

Git is installable via most systems' package managers. If you have `brew`
installed on macOS, you should already have `git`.

If you need to install `git` on Ubuntu, you can use APT:


```bash
sudo apt update
sudo apt install git
```

#### wget/cURL

`wget` and `curl` are also available from most systems' package managers. On
Ubuntu, you can use APT to install `wget`. On macOS, `curl` should already be
installed.

```bash
sudo apt update
sudo apt install wget
```

#### sratoolkit

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
    ```zsh
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

#### download_sra

!!! note
    If you downloaded this software from Zenodo, the test data SRA files are
    included in the software release zip. If you don't want to download them
    yourself, you don't need `download_sra` and can skip this step.

Make sure the `rna_clique_venv` environment created during the RNA-clique
installation is activated. Then, install dependencies for `download_sra`.

```bash
python -m pip install lxml requests
```

Clone the `download_sra` Git repository.

```bash
git clone https://github.com/actapia/download_sra
```

Add the repository root to your `PATH`.

```bash
export PATH="$PATH:$PWD/download_sra"
```

#### SPAdes

!!! note
    If you downloaded this software from Zenodo, the assembled transcriptomes
    are included in the release zip. If you don't want to assemble the
    transcriptomes yourself, you can skip downloading SPAdes.

Download the SPAdes assembler. As of this writing, 4.2.0 is the newest
version.

=== "Ubuntu"
    ```bash
    wget https://github.com/ablab/spades/releases/download/v4.2.0/SPAdes-4.2.0-Linux.tar.gz
    ```
=== "macOS"
    ```zsh
    curl -L -O https://github.com/ablab/spades/releases/download/v4.2.0/SPAdes-4.2.0-Darwin-$(uname -m).tar.gz
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

## Downloading the RNA-clique repository

Although we've already installed the RNA-clique software, a few additional files
from RNA-clique's repository will be needed for this tutorial. If you don't
already have a copy of the RNA-clique source code, clone the GitHub repository
with `git`.

<!--{{clone_command(git_branch()) | code_fence("zsh") | comment_surround}}{{empty("-->
```bash
git clone --recurse-submodules https://github.com/actapia/rna_clique
```
<!--")}}-->

We will want to keep track of the root of the RNA-clique repository in a
variable. To do this easily, change to the RNA-clique directory and then set
`RNA_CLIQUE` to `$PWD`. Finally, switch back to the parent directory. Note that
the name of the RNA-clique directory might differ if you obtained RNA-clique by
some other method than the one described above.

```bash
cd rna_clique
export RNA_CLIQUE="$PWD"
cd ..
```


## Creating a directory for our work

We will organize our files for this tutorial nicely into a single
directory. Make sure you are in the directory where you want to put the tutorial
directory, and then run the commands below to create the tutorial directory and
keep track of its path in the `TUTORIAL_DIR` environment variable.

```bash
mkdir tutorial
cd tutorial
export TUTORIAL_DIR=$PWD
```


## Obtaining sequence data

!!! note 
    If you downloaded this software from Zenodo, you already have the SRA
    files in `test_data/sra` under the root of the repository. Instead of
    completing this step, you can extract the provided data by running
    
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

We will start out with RNA-seq reads for six samples of tall fescue. These
samples are a subset of the sixteen used in the paper *RNA-clique: A method for
computing genetic distances from RNA-seq data*. The sample metadata shown
previously in the [Background](#background) section is also present in the file
[`docs/tutorials/reads2tree/tall_fescue_accs.csv`](./tall_fescue_accs.csv) in
the RNA-clique repository. We can easily download all of the RNA-seq data we
need using the `tall_fescue_accs.csv` file and the `download_sra` tool.

First, change to your `$TUTORIAL_DIR`.

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
    
    === "Ubuntu"
        ```bash
        mkdir out
        cd out
        wget "http://rna-clique-data.s3-website.us-east-2.amazonaws.com/transcripts.tar.gz"
        tar xzvf transcripts.tar.gz
        rm transcripts.tar.gz
        cd ..
        ```
    === "macOS"
        ```bash
        mkdir out
        cd out
        curl -L -O "http://rna-clique-data.s3-website.us-east-2.amazonaws.com/transcripts.tar.gz"
        tar xzvf transcripts.tar.gz 
        rm transcripts.tar.gz
        cd ..
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

Each assembled transcriptomes will be located at `transcripts.fasta` in
a directory corresponding to its sample name under the `out` directory.

## Running RNA-clique

Previous tests with this data revealed that $n = 50000$ is a good setting, so we
will use that value.

```bash
rna-clique -O "$TUTORIAL_DIR/rna_clique_out" \
           -n 50000 \
           "$TUTORIAL_DIR"/out/*
```

Verify that the `distance_matrix.h5` file was created in the output directory.

```bash
ls "$TUTORIAL_DIR/rna_clique_out/distance_matrix.h5"
```

## Results

Once RNA-clique has completed, you will likely want to [count the ideal
components](#counting-ideal-components), [see the distance
matrix](#viewing-the-distance-matrix), and possibly use the matrix as input to
downstream analyses, such as [construction of a phylogenetic
tree](#getting-a-tree) or a [PCoA plot](#getting-a-pcoa-plot), or a more direct
visualization of the matrix like a [heatmap](#getting-a-heatmap).

### Counting ideal components

The number of ideal components in an analysis is the number of genes on which
the distances are based. Since the analysis completed, we know we have at least
one ideal component, but we should hopefully have many more. An analysis based
on a number of ideal components in the single or low double digits might be too
error prone.

To check the number of ideal components, use the
[`plot_component_sizes`](../../usage.md#plot_component_sizes) program, which can
report statistics about gene matches graph components, among other things.

```bash
python -m rna_clique.plot_component_sizes --statistics \
                                          -A "$TUTORIAL_DIR/rna_clique_out"
```

You should see that we obtained around $9848$ ideal components, which is plenty
for our analysis.

### Viewing the distance matrix

To simply view the distance matrix, use the [`export_matrix`
program](../../usage.md#export_matrix).

```python
python -m rna_clique.export_matrix -O "$TUTORIAL_DIR/rna_clique_out"
```

The distance matrix should look something like this:

```text
0.0 0.0071998839240669035 0.009704850162794031 0.009340099262568933 0.009252514087212184 0.009505858709940489
0.0071998839240669035 0.0 0.009841755746060485 0.009670137702344005 0.009524426325509126 0.009712599413466836
0.009704850162794031 0.009841755746060485 0.0 0.009721949374373574 0.009786716141556896 0.009817204331389039
0.009340099262568933 0.009670137702344005 0.009721949374373574 0.0 0.009467437205501427 0.009449106568039713
0.009252514087212184 0.009524426325509126 0.009786716141556896 0.009467437205501427 0.0 0.0072939623981727415
0.009505858709940489 0.009712599413466836 0.009817204331389039 0.009449106568039713 0.0072939623981727415 0.0
```

Note that the rows and columns are not labeled. To get labels for both the rows
and columns, provide the `--format table` and `--header` options.

```python
python -m rna_clique.export_matrix --format table \
                                   --header \
                                   -O "$TUTORIAL_DIR/rna_clique_out"
```

### Downstream analyses

Although RNA-clique does not seek to integrate code for every possible
downstream analysis that could be performed on a distance matrix, and it is
assumed that many users of RNA-clique will prefer to export the matrix and use
other software for these analyses, RNA-clique nevertheless does provide a
handful of Python functions that are useful for creating trees, PCoA plots, and
heatmaps.

Since RNA-clique's downstream analysis utility functions rely on metadata about
the samples, and metadata could be expressed in a variety of formats, RNA-clique
currently does not expose the downstream analysis functions via a command-line
interface. To take advantage of the visualization functions RNA-clique provides,
one must write code that calls the visualization functions. The sections below
provide sample code that works for the data in this particular tutorial, but the
provided code might also be useful as a template for creating similar
visualizations with custom data.

#### Getting a tree

If you want a tree, you can create one using RNA-clique and Biopython. The code
below, also found in 
{{file_link("`docs/tutorials/reads2tree/make_tree.py`", "docs/tutorials/reads2tree/make_tree.py")}},
loads the
distance matrix from `distance_matrix.h5` and constructs a tree using the
neighbor-joining algorithm. The tree is also rooted at its midpoint. The tree is
saved to `nj_tree.tree`, and a visualization is saved to `nj_tree.svg` in the
`rna_clique_out` directory.

```python
--8<-- "docs/tutorials/reads2tree/make_tree.py"
```

You can run the script as follows:

```bash
python "$RNA_CLIQUE/docs/tutorials/reads2tree/make_tree.py"
```

The script creates a file called `nj_tree.svg` in the
`$TUTORIAL_DIR/rna_clique_out` directory. The result should look something like
this:

![A phylogram with six leaves representing the samples analyzed in this
tutorial. Clades containing all samples of a given genotype and only samples of
that genotype are color coded and labeled with calipers on the right side of the
figure. CTE46, CTE27, NTE, and FATG4 are in orange, blue, red, and green,
respectively.](../../images/nj_tree.svg)

#### Getting a PCoA plot

We can use the `pcoa` module to create a PCoA plot from our distance matrix. (In
turn, `pcoa` uses [`scikit-bio`](https://scikit.bio/index.html).)  To
distinguish points by genotype, we will need to use the metadata for the samples
stored at `$RNA_CLIQUE/docs/tutorials/reads2tree/tall_fescue_accs.csv`.

The code below draws a 3D and 2D PCoA plot and stores the results as SVG files
in the `rna_clique_out` directory as `pcoa_3d.svg` and `pcoa_2d.svg`,
respectively. The code can also be found at
{{file_link("`docs/tutorials/reads2tree/make_pcoa.py`", "docs/tutorials/reads2tree/make_pcoa.py")}}.

```python
--8<-- "docs/tutorials/reads2tree/make_pcoa.py"
```

The example can be run as follows:

```bash
python "$RNA_CLIQUE/docs/tutorials/reads2tree/make_pcoa.py" \
       "$TUTORIAL_DIR"/rna_clique_out
```

The two-dimensional PCoA plot should look something like this:

![A two-dimensional PCoA plot visualizing genetic distances for the six samples
used in this tutorial. Samples separate according to genotype; the two CTE27 and
CTE46 samples cluster together. Samples are color-coded by genotype. Blue points
represent CTE27, orange points CTE46, green points FATG4, and red points
NTE. The principal component axes are labeled with their relative contributions,
measured as the percentage of the sum of eigenvalues of the distance
matrix.](../../images/pcoa_2d.svg) 

The three-dimensional PCoA plot should look something like this:

![A three-dimensional PCoA plot visualizing genetic distances for the six
samples used in this tutorial. Samples separate according to genotype; the two
CTE27 and CTE46 samples cluster together. Samples are color-coded by
genotype. Blue points represent CTE27, orange points CTE46, green points FATG4,
and red points NTE. The principal component axes are labeled with their relative
contributions, measured as the percentage of the sum of eigenvalues of the
distance matrix.](../../images/pcoa_3d.svg)

#### Getting a heatmap

We can use the `draw_heatmap` function of RNA-clique to display a similarity or
distance matrix as a heatmap. The function uses the Seaborn `heatmap` function
behind the scenes, and arbitrary arguments given to `draw_heatmap` will be
passed to Seaborn.

The code below is also found in 
{{file_link("`docs/tutorials/reads2tree/make_heatmap.py`", "docs/tutorials/reads2tree/make_heatmap.py")}}.
It
draws a heatmap and saves the resulting figure in the `rna_clique_out` directory
as `distance_heatmap.svg`. 

```py
--8<-- "docs/tutorials/reads2tree/make_heatmap.py"
```

To generate a heatmap using this code, you can run the Python script as follows:

```bash
python "$RNA_CLIQUE/docs/tutorials/reads2tree/make_heatmap.py"
```

The resulting heatmap will be saved at `$TUTORIAL_DIR/distance_heatmap.svg` and
should look something like this:

![A heatmap showing distances for the six samples analyzed in this tutorial. The
heatmap is organized as a grid, and the indices shown on the left and bottom of
the heatmap indicate for each cell which pair of samples the distance shown
corresponds to. Samples are ordered and grouped by genotype on both axes. Cell
colors follow the colormap shown on the right, which maps values from $0.0095$
to $0.0075$ to colors on a gradient from dark indigo to light green. Cells
additionally show the distance values in
ten-thousandths.](../../images/distance_heatmap.svg).
