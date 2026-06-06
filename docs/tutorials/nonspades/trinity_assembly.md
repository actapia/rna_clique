# Assembling tall fescue transcriptomes with Trinity

This supplemental tutorial explains how assemble the RNA-seq reads for the six
tall fescue samples from the end-to-end ["From RNA-seq reads to a phylogenetic
tree with RNA-clique"](../reads2tree/README.md) tutorial into transcriptomes
using the [Trinity
assembler](https://github.com/trinityrnaseq/trinityrnaseq). This tutorial is a
companion to the ["Using RNA-clique with non-SPAdes data"](README.md) tutorial,
and since that tutorial provides download links for pre-assembled
transcriptomes, this assembly tutorial should be considered optional.

## Setup

### Installing Trinity

Trinity provides an [installation
guide](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Installing-Trinity),
but a different explanation of how to install Trinity via package managers is
also provided here.

!!! note
    Your system's package manager might not have the newest version of
    Trinity.

On Ubuntu, the package to install is `trinityrnaseq`. On macOS, Trinity is not
available through Homebrew's core tap but is available through
[`brewsci/bio`](https://github.com/brewsci/homebrew-bio) as `trinity`.

=== "Ubuntu"
    ```bash
	sudo apt install trinityrnaseq
	```	
=== "macOS"
    ```zsh
	brew tap brewsci/bio
	brew install trinity
	```
	
The newly installed Trinity program can be used via the `Trinity` command.

### Installing GNU Parallel (Optional)

If you have not already, consider [installing GNU
Parallel](../reads2tree/README.md#gnu-parallel-optional). GNU Parallel will make
it easier to run Trinity for multiple sets of data efficiently on multi-core
machines.

## Downloading RNA-seq reads

The process of downloading the tall fescue RNA-seq reads from the NCBI Sequence
Read Archive is described in the [Obtaining sequence
data](../reads2tree/README.md#obtaining-sequence-data) section of the end-to-end
tutorial. That section also requires that some previous sections of the tutorial
be completed. Specifically, it is expected that the
[Setup](../reads2tree/README.md#setup) and [Creating a directory for our
work](../reads2tree/README.md#creating-a-directory-for-our-work) sections have
been completed.

It will be assumed that the `SRR*.fastq` files downloaded in the end-to-end
tutorial are in the `$TUTORIAL_DIR` for the remainder of this tutorial.

## Assembling transcriptomes

First, change to the `$TUTORIAL_DIR`.

```bash
cd "$TUTORIAL_DIR"
```

We will run Trinity with multiple CPU cores and a memory limit. What specific
settings you should use here depends on what hardware you have available on your
machine. 

### CPU usage

The product of the number of Parallel jobs and the number of CPUs used
by Trinity should not exceed the number of logical cores you have on your
system, and you should prefer to use more Parallel jobs rather than more CPUs
with Trinity. If you're not using Parallel, then the number of CPUs used by
Trinity must not exceed the number of logical cores you have.

### Memory usage

The product of the number of Parallel jobs and memory limit should not exceed
the amount of RAM you have in your system. If you are not using Parallel, then
the memory limit must not exceed the amount of RAM you have in your system.

The examples below are written for a machine with 32 cores and over 300 GB of
RAM. You will need to change the parameters to fit your own system.

=== "With Parallel"
```bash
parallel -j 6 Trinity --seqType fq --max_memory 48G --single {} --CPU 2 \
         --output trinity_{/.} ::: *.fastq
```
=== "Without Parallel"
```bash
for f in *.fastq; do
    s="${f%%.*}";
	Trinity --seqType fq --max_memory 256G --single "$f" --CPU 32 \
	        --output "trinity_$s";
done
```
