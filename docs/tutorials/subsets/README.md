# Quickly computing subsets of existing analyses

This tutorial explains how (and why) to use the
[`make_subset.py`](usage.md#make_subsetpy) script to create subsets of existing
RNA-clique analyses. This guide should be useful to anyone who wishes to analyze
a subset of the samples analyzed in a previously completed RNA-clique analysis.

This tutorial will assume that the user is starting from the end of the
end-to-end ["From RNA-seq reads to a phylogenetic tree with
RNA-clique"](../reads2tree/README.md) tutorial; this will give us a completed
analysis to work from.

## Check your environment

Before we begin, check that your `RNA_CLIQUE` and `TUTORIAL_DIR` environment
variables are set and pointing to the RNA-clique root directory and the working
directory for the tutorial, respectively.

```bash
echo "$RNA_CLIQUE"
echo "$TUTORIAL_DIR"
```

If either of these commands prints an empty line or prints a path to a directory
that does not exist, you may need to reset the environment variables to the
values described in the ["Creating a directory for our
work"](../reads2tree/README.md#creating-a-directory-for-our-work) section of the
end-to-end tutorial.

Also, check that you are in the `rna-clique` Conda environment. You should see
`(rna-clique)` at the beginning of your prompt. If not, try activating the
environment with

```bash
conda activate rna-clique
```

## Goal and rationale

For reference, the genotypes and endophyte statuses of the samples analyzed in
the "From RNA-seq reads to a phylogenetic tree with RNA-clique" tutorial are
shown below:

| SRA Accession | Genotype | Endophyte |
|---------------|----------|-----------|
| SRR2321388    | CTE46    | infected  |
| SRR2321385    | CTE46    | minus     |
| SRR8003761    | CTE27    | infected  |
| SRR8003762    | CTE27    | minus     |
| SRR7990321    | FATG4    | infected  |
| SRR8003736    | NTE      | infected  |

As explained in the end-to-end tutorial, some of the samples possess an
endosymbiotic fungus, *Epichloë coenophiala*. The presence or absence of this
fungus each of the samples is shown in the "Endophyte" column.

Suppose we wanted to look only at the infected samples. (We don't actually
expect to see anything interesting among the infected samples specifically, but
this subset selection is good for a tutorial and will be useful for a later
tutorial, "[Exporting and searching ideal
components](../export_and_search/README.md)"). Of course, we could just take the
rows and columns corresponding to the infected samples from of the distance
matrix from our full analysis. In general, this approach is not the best option
because we can often obtain more ideal components with smaller subsets of
samples, and simply taking the distances directly from the existing distance
matrix forces us to use the ideal components from the full (supserset)
analysis. With more ideal components, we generally expect more precise
results. Additionally, very distant pairs of samples might also makes distances
between very closely related pairs of samples less accurate. This kind of
problem is most pronounced when the full analysis consists of two or more tight
clusters that are distant from each other, and the subset analysis considers
only one of these clusters. For example, if the full analysis considers two
species simultaneously, we might expect two distant species clusters, and simply
taking the distances for one species out of the full analysis's distance matrix
might give relatively imprecise distances.

Actually, for this particular case, obtaining the infected-only distance matrix
directly from the full distance matrix might be a reasonable option. We have
plenty of ideal components, and since the minus samples do not cluster together,
they are not so distant from the infected genotypes that the distances among the
infected genotype samples are obscured. For the sake of example, however, we'll
assume that we need to recompute the distances for the set of infected
samples. Although the recomputation may not be necessary, it would be
interesting to see how many more ideal components we obtain when analyzing only
the infected samples.


A naive way of recomputing a distance matrix for just the infected samples would
be to re-run the command in the ["Running
RNA-clique"](../reads2tree/README.md#running-rna-clique) section while replacing
the last argument with
`"$TUTORIAL_DIR"/out/{SRR2321388,SRR2321385,SRR8003761,SRR8003762}`. Although
this approach is straightforward, it be inefficient because we already did most
of the work needed to recompute this distance matrix when we ran the original
subset analysis. If we were to re-run RNA-clique with just the infefcted
samples, RNA-clique would need to first select the top $n$ genes from the
infected samples, then obtain gene matches tables for each pair of infected
samples, and these steps would account for the majority of the running time of
RNA-clique. 

In the superset analysis, however, RNA-clique also selected the top $n$ genes
from the infected samples (along with all of the other samples) and computed the
gene matches tables for every pair of infected samples (and every other pair of
samples). The idea behind the "fast subsetting" feature provided by the
`make_subset.py` script is to reuse the existing gene matches tables and top
genes files for the new subset analysis. From the gene matches tables for just
the pairs of infected samples, `make_subset.py` can build a new gene matches
graph, filter the gene matches tables, and compute new pairwise distances, and
these steps take little time compared to the full RNA-clique method.

## Creating a subset analysis

We will use the `make_subset.py` script to create a subset of the previous
analysis. `make_subset.py` allows the user to select a subset by specifying
criteria for samples to include or exclude from the subset. `make_subset.py`
allows exclusion criteria to be specified via a list of samples to exclude,
which may be given as arguments in the command line invocation of
`make_subset.py` or in a file provided to `make_subset.py`. Inclusion criteria
can be specified via a list of samples to include, likewise provided as either
arguments in the command line invocation or in a specified file, but inclusion
criteria can also be specified using a regular expression that matches only
samples to include. For this example, we'll specify the infected samples
directly on the command line.

In the `make_subset.py` command, we must also tell the script which analysis we
are subsetting, specified via a path to the configuration file for the analysis,
and the output directory in which the subset analysis files will reside.

```bash
python make_subset.py -I "$TUTORIAL_DIR"/rna_clique_out/config.yaml \
                      -O "$TUTORIAL_DIR"/infected_subset_out \
                      -y SRR2321388 SRR8003761 SRR7990321 SRR8003736
```

You should find that the above command completes *much* more quickly than a full
RNA-clique run would, but if you check the subset output directory at
`"$TUTORIAL_DIR"/infected_subset_out`, you'll see that we don't have a distance
matrix yet. To get the distance matrix, we also need to run
`filtered_distance.py` on the subset analysis.

```bash
python filtered_distance.py -O "$TUTORIAL_DIR"/infected_subset_out
```

If you check `"$TUTORIAL_DIR"/infected_subset_out` again, you should now see a
`distance_matrix.h5` file.

We can also check the number of ideal components.

```bash
python plot_component_sizes.py --statistics \
                               -A "$TUTORIAL_DIR"/infected_subset_out
```

You should see that we have more ideal components (around $12676$) with just the
four infected samples.
