# Exporting and searching ideal components

The ideal components constructed by RNA-clique can be interpreted as sets of
orthologous genes identified among all of the analyzed samples. RNA-clique uses
alignments among these identified orthologs as the basis for its distance
computation without specifically selecting and writing the ortholog sequences to
the disk, but since it is sometimes useful to be able to manually inspect or
search the detected orthologs, RNA-clique can optionally export and search
orthologous genes in ideal components via the `export_and_search.py`
script. This tutorial describes how to use the `export_and_search.py` script to
export orthologous gene sequences from ideal components and search those genes
for a sequence of interest.

The first part of this tutorial assumes that the user has completed the
end-to-end ["From RNA-seq reads to a phylogenetic tree with
RNA-clique"](../reads2tree/README.md) tutorial, and the second part assumes that
the user has also completed the ["Quickly computing subsets of existing
analyses"](../subsets/README.md) tutorial. These earlier tutorials give us some
data to export and search.

## Background

The metadata for the six samples analyzed in the end-to-end tutorial is
reproduced below.

| SRA Accession | Genotype | Endophyte |
|---------------|----------|-----------|
| SRR2321388    | CTE46    | infected  |
| SRR2321385    | CTE46    | minus     |
| SRR8003761    | CTE27    | infected  |
| SRR8003762    | CTE27    | minus     |
| SRR7990321    | FATG4    | infected  |
| SRR8003736    | NTE      | infected  |

Recall that some of the samples have an endophyte&mdash;the fungus *Epichloë
coenophiala*&mdash;such samples are those labeled as "infected" in the
endophyte column of the sample metadata. In this tutorial, I will also describe
such samples as "E+" ("E plus") and will refer to samples without endophytes as
"E&minus;" ("E minus"). In this exercise, we will export the orthologs from the
full six-sample analysis and search for sequences from the genome of *Epichloë
coenophiala* within the ideal components. We will then repeat the process for
the subset of four infected samples analyzed in the subsetting tutorial.

We know that the sequence reads from which we assembled the E+ samples also
contain RNA sequences from the endophyte, and, in fact, some of these reads were
assembled into transcripts by rnaSPAdes. Since the E&minus; samples lack
endophytes, there should be neither endophyte RNA-seq reads nor assembled
transcripts for those samples.

Endophyte transcripts in the assembled transcriptomes could potentially
introduce error into the distance calculations if alignments between endophyte
transcripts were used in the distance calculation. For such transcripts to be
used, their genes would need to appear in some ideal component in the gene
matches graph. In the case of the full six-sample analysis, we expect that no
such endophyte genes will be in ideal components because endophyte transcripts
should not have identifiable orthologs in the E&minus; transcriptomes. For the
subset of four infected samples, however, the presence of endophyte genes in all
samples means there is a chance that some of the genes will end up in ideal
components (i.e., have orthologs in all samples). In this tutorial, we will
check this intuition using RNA-clique's export and searching features.

## Check your environment variables

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

## Download Epichloë genome assemblies

We will download two genome assemblies of *Epichloë coenophiala* to use as query
sequences for searching the exported ideal components in both the original and
the infected subset RNA-clique analyses. First, create a directory for the
genomes

```bash
mkdir "$TUTORIAL_DIR"/ec_genomes
```

To download the genomes, run

```bash
wget -P "$TUTORIAL_DIR"/ec_genomes/ \
        https://cs.uky.edu/~acta225/rna_clique/ec_genomes.tar.xz
```

Then, extract the genomes with `tar`.

```bash
tar xJvf "$TUTORIAL_DIR"/ec_genomes/ec_genomes.tar.xz \
         -C "$TUTORIAL_DIR"/ec_genomes
```

Verify that the genomes have been extracted to `$TUTORIAL_DIR/ec_genomes`.

```bash
ls "$TUTORIAL_DIR"/ec_genomes/
```

You should see three files: `ec_genomes.tar.xz`,
`e19_scaffolds-JAFEMN000000000.fasta`, and
`e4305_Mas339_20200623_Ref_Scaffolds_CLS.fasta`.

## Export and search ideal components in the original analysis

To export and search the orthologous genes in ideal components at once, run the
`export_and_search.py` script. We'll use a low e-value threshold, `1e-99`, to
ensure that we only get very close matches.

```bash
python export_and_search.py -C "$TUTORIAL_DIR"/rna_clique_out/config.yaml \
                            -Q "$TUTORIAL_DIR"/ec_genomes/*.fasta \
                            -X "$TUTORIAL_DIR"/full_ec_search_out \
                            -e 1e-99
```

The results should be written to
`$TUTORIAL_DIR/full_ec_search_out/rna_clique_out`. (If we had provided
`export_and_search.py` multiple configuration files, there would be multiple
directories under `$TUTORIAL_DIR/full_ec_search_out/`.) 

## Examining and interpreting results

### Exported ideal components

The exported ideal components should be under the
`$TUTORIAL_DIR/full_ec_search_out/rna_clique_out/export` directory. If you list
the contents of the `export` directory, you should see many FASTA files
corresponding to the ideal components identified by RNA-clique. If you open one
of these file, you will see that the FASTA headers differ from those in the
original input transcriptome files. For example, a sequence might have the FASTA
header `>-NODE_1_length_14835_cov_23.346544_g0_i0:SRR7990321`. In this modified
header, the leading `-` denotes that the sequence is the reverse complement of
the original sequence in the input transcriptome. The sequence was reoriented to
ensure that all sequences in the same file are in the same orientation. The
trailing `:SRR7990321` indicates that the original sequence comes from the
transcriptome of sample `SRR7990321`. Everything between the leading `-` and
trailing `:` is the original FASTA sequence header from the input transcriptome.

Other FASTA headers may lack the leading `-`. In that case, the sequence appears
exactly as it did in the input transcriptome; it has not be reoriented.

Under the `export` directory is also a file called `all_ideal.fasta`, which
contains all sequences from the exported ideal components. `all_ideal.fasta` is
akin to a concatenation of all of the `ideal_component` FASTA files, but
sequence headers in the `all_ideal.fasta` are further modified to indicate from
which ideal components they come. The ideal component to which a sequence
belongs is indicated by a suffix appended to the sequence header. For example,
the header
`>NODE_1_length_15383_cov_32.255511_g0_i0:SRR2321385:ideal_component_0` comes
from `ideal_component_0`.

### Search results

The search results are
organized by query file, so we should have directories named
`search_e19_scaffolds-JAFEMN000000000` and
`search_e4305_Mas339_20200623_Ref_Scaffolds_CLS` under
`$TUTORIAL_DIR/full_ec_search_out/rna_clique_out`. Since we did not perform an
"extended search" to find higher e-value matches in ideal components where there
were already low e-value matches, we should only see three files in each of
these two directories&mdash;`queries.sam`, `stats`, and `subjects.fasta`. (For
an explanation of the extended search and the files produced in such a search,
see the Command-line usage guide entry for
[`search_deail_components.py`](../../usage.md#search_ideal_componentspy).)

#### stats

The `stats` file provides statistics for the search results in JSON format. Each
`stats` file should have three JSON keys&mdash;`hits`, `seqs`, and
`components`. `hits` is the number of BLAST HSPs found in the search. `seqs` in
the number of distinct transcripts matches in the search, and `components` is
the number of ideal components in which at least one transcript matched.

For `search_e19_scaffolds-JAFEMN000000000` directory, the `stats` file should
look something like this:

```json
{"hits": 124, "seqs": 43, "components": 15}
```

For the `search_e4305_Mas339_20200623_Ref_Scaffolds_CLS` directory, the `stats`
file should look something like this:

```json
{"hits": 87, "seqs": 44, "components": 16}
```

Surprisingly, in both cases, several components contained matching sequences,
even though some of the input samples possessed no endophyte. It's possible,
however, that the matching transcripts were all from the E+ samples, and that
the matches are simply misassembled transcripts that happened to incorporate
reads from the endophyte. Alternatively, it could be that the query genome
assemblies inadvertently incorporated reads from the host plants. We will
investigate these possibilities later in this tutorial.

#### subjects.fasta

The `subjects.fasta` file contains the sequences of all transcripts from the
`all_ideal.fasta` file that matched a query sequence in the BLAST search. The
number of sequences in this file should match the value associated with the
`seqs` key in the `stats` file.

Since the sequences in `subjects.fasta` are taken from `all_ideal.fasta`, each
sequence header indicates the sample and ideal component in which the transcript
is found.

#### queries.sam

`queries.sam` contains the BLAST alignments of the query sequences to the ideal
component transcript sequences (in `subjects.fasta`) in [Sequence Alignment
Map](https://samtools.github.io/hts-specs/SAMv1.pdf) (SAM) format. The number of
alignments in this file should match the value associated with the `hits` key in
the `stats` file.

### Checking represented samples in matching components

As mentioned in the [`stats`](#stats) subsection, the endophyte sequences found
among ideal components could possibly be explained by transcriptome assembly
errors. Specifically, it is possible that RNA-seq reads from the endophyte were
erroneously incorporated into plant transcripts, creating spurious *in silico*
hybrid transcripts. If that were the case, then we would expect all hits to
occur in E+ sample where such endophyte reads would be present.

We use the sequence headers in `subjects.fasta` to check whether matches all
come from E+ samples.

```bash
grep "$TUTORIAL_DIR"/full_ec_search_out/*/search_*/subjects.fasta \
     -e '>' | grep '[^:]*:[^:]*$' -o | sort -u | sort -t: -k2
```

The command above gets the sample name and ideal component from every FASTA
header in `subjects.fasta`. `sort` is then used twice to get the unique
sample-ideal component pairs and then sort them by ideal component.

You should see something like this:

```text
SRR7990321:ideal_component_1127
SRR8003761:ideal_component_1331
SRR8003761:ideal_component_1459
SRR2321388:ideal_component_2366
SRR2321388:ideal_component_3289
SRR8003761:ideal_component_3939
SRR2321388:ideal_component_4110
SRR8003761:ideal_component_511
SRR8003761:ideal_component_6114
SRR2321385:ideal_component_6249
SRR2321388:ideal_component_6249
SRR7990321:ideal_component_6249
SRR8003736:ideal_component_6249
SRR8003761:ideal_component_6249
SRR8003762:ideal_component_6249
SRR8003761:ideal_component_667
SRR2321385:ideal_component_6777
SRR2321388:ideal_component_6777
SRR7990321:ideal_component_6777
SRR8003736:ideal_component_6777
SRR8003761:ideal_component_6777
SRR8003762:ideal_component_6777
SRR7990321:ideal_component_7010
SRR8003761:ideal_component_7855
SRR2321388:ideal_component_793
SRR7990321:ideal_component_8416
```

Notice that transcripts from SRR2321385 in ideal components $6249$ and $6777$
matched the endophyte assemblies. We can conclude that the hits are not merely a
result of assembly errors incorporating endophyte reads into plant sequences.

### Visualizing alignments produced by export\_and\_search.py

Since `export_and_search.py` produces alignments in SAM format, we can use a
standard genome browser to view the alignments. In this example, I will show the
alignments using [Integrative Genomics Viewer](https://igv.org/doc/desktop/)
(IGV).

Before visualizing the alignments, make sure that you have copies of the
`subjects.fasta` and `queries.sam` files on the computer that will be running
IGV, and make sure you can find these files on your system.

#### Downloading, installing, and launching IGV

IGV can be installed from https://igv.org/doc/desktop/#DownloadPage/ . Select
the download most appropriate for your system. In most cases, you will want to
install a version with Java included. (Otherwise, IGV will try to use the system
Java installation, which may not be compatible.)

On Windows, run and complete the installer. You should be able to find IGV in
your start menu.

On macOS, open the `.zip` file you downloaded and move the app within to your
Applications folder, then launch IGV from there.

On Linux, unzip the downloaded archive and run `igv.sh`.

#### Loading the subjects.fasta file

Before we can load the alignments, we need to load the reference sequences into
IGV. IGV considers any reference sequences to be "genomes," but it will work
with our transcripts as well.

In the menu bar, under **Genomes**, select "Load Genome from File...". A file
selection dialog will appear; select the `subjects.fasta` file from the
`e19_scaffolds-JAFEMN000000000` search.

#### Loading the queries.sam file

Now that the references (subjects) have been loaded, we can open alignments to
display as features along the subjects.

In the menu, under **File**, select "Load from File..." and choose the
`queries.sam` file from the file selection dialog that appears.

You will most likely be prompted about a missing index file for the
`queries.sam` alignments. Click "Go" to create the index file automatically.

#### Viewing the alignments

By default, IGV is set to try to view "All chromosomes," i.e., all subject
sequences, simultaneously. In this mode, no alignments will appear, so we need
to select a "chromosome" (subject sequence) from the drop down menu at the top
of the window. Let's select the first transcript from ideal component 6249. You
should see alignments like the ones below:

![Two segments of an Epichloë coenophiala genome shown aligned to a transcript
from an assembled transcriptome for tall fescue. The two alignments are visually
very similar, though they differ in some gap and insertion locations. Both
alignments span the majority of the transcript.](igv_snapshot.svg)

You should see the alignments span most of the transcript, and if you look at
the other transcripts from the same ideal component, you should see the same
thing for the other alignments. This is true even for the E&minus; samples, such
as SRR2321385, and the alignments clearly span much more than the length of a
single read.

#### Searching the Nucleotide database for matches transcripts

One other possible explanation for the alignments we see is that the genome
assemblies could inadvertently have included some sequences from the host
plant. We can check if this is true by BLASTing the transcripts in the matched
exported ideal components against the NCBI non-redundant Nucleotide
database. From such a BLAST search, we could see in which kinds of organisms we
get hits. We will also learn more about the identity of the matching sequences,
which might help explain the matches to the E&minus; samples. 

First, let's concatenate the exported transcripts from the ideal components
where all samples had matching transcripts. For this run, we want ideal
components $6249$ $6777$, but beck the results you got at the end of the
["Checking represented samples in matching
component"](#checking-represented-samples-in-matching-components) section to
make sure you are using the right ideal component IDs here.

```bash
cat "$TUTORIAL_DIR"/full_ec_search_out/*/export/ideal_component_6249.fasta \
    "$TUTORIAL_DIR"/full_ec_search_out/*/export/ideal_component_6777.fasta \
	> matching_components.fasta
```


Now, we will search the created `matching_components.fasta` against the NCBI
Nucleotide database. (You may prefer to do this step in the web BLAST interface
to avoid rate limits.)

```bash
blastn -query matching_components.fasta -evalue 1e-99 -remote -db nr -outfmt 6 \
       -out remote_results
```

If you view the subject sequence IDs in the resulting `remote_results` file and
search for them within the
[Nucleotide](https://www.ncbi.nlm.nih.gov/nucleotide/) database, you should find
that 

Actually, further investigation of the sequences we found reveals that they code
for highly conserved proteins, and we are getting matches to both a plant and
the fungal homolog of the gene. This explains why we get matches for the
endophyte sequence even in E&minus; samples.

## Export and search ideal components in the subset analysis

Now that we've looked at the results for the full analysis, we will see how our
results differ for the subset analysis that included only E&minus; samples. Run
`export_and_search.py` again but provide it with the configuration file for the
subset analysis and a different output directory.

```bash
python export_and_search.py -C "$TUTORIAL_DIR"/infected_subset_out/config.yaml \
                            -Q "$TUTORIAL_DIR"/ec_genomes/*.fasta \
                            -X "$TUTORIAL_DIR"/subset_ec_search_out \
                            -e 1e-99
```

## Examining results for the subset analysis

The output export and search files should be located at
`$TUTORIAL_DIR/subset_ec_search_out/infected_subset_out`. Check the `stats`
files under `search_e19_scaffolds-JAFEMN000000000` and
`search_e4305_Mas339_20200623_Ref_Scaffolds_CLS`. You should see that the number
of components in which matching transcripts have been found has increased
dramatically compared to the results for the original set of samples, from
around `15` to over `500` for both genomes. Although the subset only had around
120% as many ideal components as the original set, the subset has at least 3200%
as many ideal components matching the endophyte genomes. This suggests that
RNA-clique filters out some of the endophyte transcripts when the full analysis
with E- samples is performed, but RNA-clique cannot filter such transcripts when
all the samples are E+.
