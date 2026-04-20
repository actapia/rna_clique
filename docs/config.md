# Configuration files

In addition to a basic [command-line interface](usage.md), most of the scripts
that are part of RNA-clique can be controlled via [YAML](https://yaml.org/)
configuration files. Configuration files can be useful for making analyses more
reproducible or for performing multiple analyses with similar parameters. Many
RNA-clique programs can also *output* configuration files to store the 
parameters of the analysis performed in a machine-readable format that can be
readily used to repeat the analysis or as a template for further analyses with
similar parameters.

The document first provides a [basic template](#template) configuration file
with default or empty values for all settings. This document then describes how
configuration settings may be [derived automatically](#derived-settings) from
others and how configuration files interact with the [command-line
interface](#interaction-with-command-line-arguments). This document discusses
[optional values](#optional-values) and how [file paths](#file-paths) are
interpreted.  Finally, this document provides a reference for the
[settings](#settings) that can be provided in RNA-clique configuration files and
explains how to [troubleshoot](#troubleshooting-configuration) RNA-clique's
configuration.

## Template

```yaml
# Version of the configuration schema used.
config_version: 0.0.1
# Name to assign to the analysis.
title:
# Directories containing the transcript FASTA files.
input_dirs:
# Output directory containing top n genes by coverage.
top_genes_dir:
# Output directory containing gene matches tables.
tables_dir:
# Intermediate directory containing BLAST DB caches.
cache_dir:
# Output directory root.
output_dir:
# Output gene matches graph.
graph:
# Number of top genes by k-mer coverate to select.
top_genes:
# Name of transcripts files in input directories.
transcripts_name: transcripts.fasta
# Threshold for counting a match between two genes.
top_matches: 1
# e-value threshold to use for BLASTn searches.
evalue: 1e-99
# Keep all matches between genes in the case of ties.
keep_all: true
# Number of parallel jobs to use.
jobs: 31
# Python regex to use for parsing transcript IDs.
transcript_id_regex: ^.*cov_([0-9]+(?:\.[0-9]+))_g([0-9]+)_i([0-9]+)
# Mapping from paths to sample names.
path_to_sample:
# Output distance matrix location.
matrix:
# When the last analysis associated with this config file finished.
finished:
# Version of RNA-clique used to create this analysis.
version: v0.2.0
# Path to analysis of which this is a subset.
subset_of:
```

## Derived settings

In some cases, RNA-clique will derive good default values for some settings
automatically from the values of other settings. Those default values can be
used if the user does not specify a value for the setting explicitly. RNA-clique
produces such default settings to simplify RNA-clique commands for common use
cases by allowing users to provide fewer command-line arguments. Although
RNA-clique allows command-line invocations simplified in this way, it retains
the flexibility gained from having many command-line arguments by still allowing
potentially derived settings to be specified explicitly.

For example, RNA-clique will automatically derive values for `top_genes_dir`,
`cache_dir`, `tables_dir`, `graph`, `matrix`, and `output_config` from
`output_dir`, which means only a root `output_dir` value must be provided. Any
of the derived default values for these settings can be overridden, however,
when the setting is defined explicitly by the user.

How settings can be derived from others depends on the script being used, though
many scripts share similar rules. These rules are documented for each script in
the [Command-line usage guide](usage.md). When a setting can be derived from
another for some script, the rule is described in the "default value" column of
the row for the derived setting in the table under the script's "Options"
subsection. In most cases, the rule is described in a kind of succinct
pseudocode, but when the rule cannot be expressed so succinctly, the table
simply documents the settings on which the derived setting depends. Rarely,
there may be multiple possible rules for deriving a value. In such cases, the
table simply lists the default value as "dynamic," and the rules are explained
in more detail beneath the table.

## Interaction with command-line arguments

RNA-clique can be controlled via command-line arguments and configuration files
simultaneously by providing both command-line arguments and a configuration file
via the `--input-config/-C` option. Using settings from both sources can make
RNA-clique more convenient to use; for example, it makes it possible to put
common settings in a configuration file and make further adjustments to settings
without requiring that the user open an editor to change the configuration file
each time the user wants to try a different setting. Unfortunately, accepting
settings from multiple sources inevitably raises questions about how those
sources should interact, and the problem is exacerbated somewhat when certain
settings can be derived automatically from other settings, as they are in many
of RNA-clique's programs. This section briefly discusses rules for using
command-line arguments and configuration files simultaneously. Of course, a user
can always opt to use only one source at a time should they prefer to avoid the
complexity of the interactions described here.

Currently, RNA-clique initially processes the arguments provided directly on the
command line. Where possible, any rules for deriving settings not provided
explicitly are applied until no further derivations are possible.

After the command-line arguments have been processed, RNA-clique loads the
settings provided in the configuration file (if specified). Settings specified
explicitly in the configuration file *do not* override any settings specified
explicitly at the command line or derived from the command-line settings. After
the explicit configuration file settings have been loaded, the rules for
deriving settings not provided explicitly are again applied to derive missing
settings. 

<!-- At this stage, the derived settings from the previous step might be overwritten, -->
<!-- depending on what settings are available and the priorities of the rules used -->
<!-- for deriving the settings. In most cases in RNA-clique, however, no such -->
<!-- overwriting will happen. -->

Since the rules above are slightly complex, it is recommended that users
manually using configuration files and command-line arguments simultaneously
take advantage of the [`--show-config` option](#troubleshooting-configuration)
to verify that RNA-clique is using the expected settings. 

## Optional values

Configuration files provided to RNA-clique programs need not provide values for
all possible settings. Values missing from the configuration can be provided
[via command-line arguments](#interaction-with-command-line-arguments) or
[derived from the provided values automatically](#derived-settings), and an
RNA-clique program will not fail if settings it does not need are missing. The
[Command line usage guide](usage.md) documents which settings must be present
(from any source) in the "Required" column of each script's "Options" table.

## File paths

When a setting is assigned a relative file path, the file path is interpreted as
relative to the working directory of the *running script* rather than the
directory of the configuration file or the analysis root. This behavior ensures
that relative paths are treated the same whether they are provided via a
configuration file or via a command-line argument, but it can cause some
confusion when reusing a config file from a different directory. Future versions
of the RNA-clique configuration schema may introduce a `working_directory`
setting to document where the analysis was executed. 

<!-- RNA-clique may require that -->
<!-- the `working_directory` be the same as the working directory from which the -->
<!-- script is running, or RNA-clique might reinterpret paths provided in the -->
<!-- configuration file based on the `working_directory` setting. -->

## Settings

| Setting                                                | Python type               | YAML type                     | Description                                                       |
|:-------------------------------------------------------|:--------------------------|:------------------------------|:------------------------------------------------------------------|
| [`config_version`](config.md#config-version)           | `str`                     | Scalar                        | Version of the configuration schema used.                         |
| `title`                                                | `str`                     | Scalar                        | Name to assign to the analysis.                                   |
| [`input_dirs`](config.md#input-dirs)                   | `list[pathlib.Path]`      | Sequence of Scalar            | Directories containing the transcript FASTA files.                |
| [`top_genes_dir`](config.md#top-genes-dir)             | `pathlib.Path`            | Scalar                        | Directory containing top n genes by coverage.                     |
| [`tables_dir`](config.md#tables-dir)                   | `pathlib.Path`            | Scalar                        | Directory containing gene matches tables.                         |
| [`cache_dir`](config.md#cache-dir)                     | `pathlib.Path`            | Scalar                        | Directory containing BLAST DB caches.                             |
| [`output_dir`](config.md#output-dir)                   | `pathlib.Path`            | Scalar                        | RNA-clique analysis output root directory.                        |
| [`graph`](config.md#graph)                             | `pathlib.Path`            | Scalar                        | Gene matches graph.                                               |
| [`top_genes`](config.md#top-genes)                     | `int`                     | Scalar                        | Number of top genes by k-mer coverate to select.                  |
| [`transcripts_name`](config.md#transcripts-name)       | `str`                     | Scalar                        | Name of transcripts files in input directories.                   |
| [`top_matches`](config.md#top-matches)                 | `int`                     | Scalar                        | Threshold for counting a match between two genes.                 |
| `evalue`                                               | `float`                   | Scalar                        | e-value threshold to use for BLASTn searches.                     |
| [`keep_all`](config.md#keep-all)                       | `bool`                    | Scalar                        | Keep all matches between genes in the case of ties.               |
| `jobs`                                                 | `int`                     | Scalar                        | Number of parallel jobs to use.                                   |
| [`transcript_id_regex`](config.md#transcript-id-regex) | `re.Pattern`              | Scalar                        | Python regex to use for parsing transcript IDs.                   |
| [`path_to_sample`](config.md#path-to-sample)           | `dict[pathlib.Path, str]` | Mapping from Scalar to Scalar | Mapping from paths to sample names.                               |
| [`matrix`](config.md#matrix)                           | `pathlib.Path`            | Scalar                        | Output distance matrix location.                                  |
| `finished`                                             | `datetime.datetime`       | Scalar                        | When the last analysis associated with this config file finished. |
| `version`                                              | `str`                     | Scalar                        | Version of RNA-clique used to create this analysis.               |
| [`subset_of`](config.md#subset-of)                     | `pathlib.Path`            | Scalar                        | Path to analysis of which this is a subset.                       |

### config\_version

RNA-clique configuration files can state what version of the configuration
schema they use. This setting is designed to provide a measure of
future-proofing in case eventual changes to the schema require explicit
detection of the schema version to provide backwards compatibility. Currently,
the schema version is `0.0.2`, but no programs in RNA-clique use the schema
version.

### input\_dirs

The `input_dirs` setting should be a list (YAML Sequence) of file paths to
[directories containing the transcriptome files](formats.md#transcriptomes).

### top\_genes\_dir

The `top_genes_dir` setting should be a path to a 
[directory containing the top genes](formats.md#top-genes) file.

### tables\_dir

The `tables_dir` setting should be a path to a [directory containing the gene
matches tables](formats.md#gene-matches-tables).

### cache\_dir

The `cache_dir` is a directory in which RNA-clique will store a `simple_blast`
[BLAST database
cache](https://github.com/actapia/simple_blast/blob/main/README.md#db-caches)
containing BLAST databases for the top genes files. `cache_dir` is an
intermediate output directory not intended to be used by the RNA-clique user.

### output\_dir

`output_dir`, when provided, is used as the default "root" for the RNA-clique
analysis. Paths for for `top_genes_dir`, `cache_dir`, `tables_dir`, `graph`,
`matrix`, and `output_config` will usually be [automatically
placed](#derived-settings) under the `output_dir` when those settings are
not provided explicitly. (Check the relevant "Options" section for the script in
the [Command-line usage guide](usage.md) to verify this behavior for the script
being used.) 

### graph

The `graph` setting should be a path to the [gene matches
graph](formats.md#gene-matches-graph). 

### top\_genes

`top_genes` is the number of top genes to select by $k$-mer coverage when
creating the [top genes files](formats.md#top-genes). In the original RNA-clique
paper, `top_genes` is referred to as parameter $n$. 

$k$-mer coverage quantifies the amount of sequence data that contributes to an
assembled transcript. Some assemblers, including
[rnaSPAdes](https://github.com/ablab/spades), report $k$-mer coverage values for
each assembled transcript. When an assembler can report multiple isoforms of a
single transcript, transcripts can be grouped into "isotig sets" representing
isoforms of the same gene. (rnaSPAdes does this grouping automatically, but a
separate program could also be used to organize transcripts into isotig sets.)
Although the assembler only assigns $k$-mer coverage values to individual
transcripts, RNA-clique takes the $k$-mer coverage of a gene to be the maximum
$k$-mer coverage among all transcripts that belong to the same gene (isotig
set).

When RNA-clique selects the top $n$ genes by $k$-mer coverage, it is intended to
select those genes best supported by the RNA-seq data. This step is intended to
both reduce errors caused by poorly assembled transcripts and reduce the amount
of time needed to perform the remainder of the analysis. What setting of $n$ is
best depends on the data, but $n = 50000$ worked well for the analyses performed
in the original RNA-clique.

### transcripts\_name

RNA-clique expects all input [transcriptome](formats.md#transcriptomes) FASTA
files to have the same filename. By default, RNA-clique expects them all to be
named `transcripts.fasta` since this is the default output filename for
[rnaSPAdes](https://github.com/ablab/spades), but this name can be customized
using the `transcripts_name` setting.

### top\_matches

When comparing sample $A$ and sample $V$, we BLAST $A$ against $B$ and $B$
against $A$. Ordinarily, we keep a pair of genes $g$ (from $a$) and $h$ (from
$b$) when merging the results from the two directions if and only if $g$ is
among the best matches for $h$ in $A$, *and* $h$ is among the best matches for
$g$ in $B$, according to bitscore. (We allow ties, so $h$ may not be the *only*
best match for $g$ in $B$, and, likewise, $g$ may not be the *only* best match
for $h$ in $A$.) This behavior corresponds to a parameter setting of
`top_matches = 1` ($N = 1$ in the original RNA-clique paper) because we are
consider only the matches with top $N = 1$ bitscore in both directions.

We could alternatively set $N$ to some value greater than 1. In that case, when
we merge the two directions, we could keep a pair of genes $g$ and $h$ if and
only if $g$ is among the matches with top $N$ bitscore for $h$ in $A$, and $h$
is among the matches with top $N$ bitscore for $g$ in $B$.

As of this writing, values for `top_matches` greater than $1$ are mostly
untested, and it is recommended that this parameter simply be set to $1$ in
practice.

### transcript\_id\_regex

The `transcript_id_regex` setting should be a Python regular expression that can
be used to parse the FASTA sequence header lines of the transcript FASTA
files. The following capture groups are expected and can be identified by
position or by
[name](https://docs.python.org/3/library/re.html#regular-expression-syntax).

|   Position | Name     | Description                                            |
|-----------:|:---------|:-------------------------------------------------------|
|          1 | coverage | $k$-mer coverage, expressed as a floating-point number |
|          2 | gene     | Gene ID, a non-negative integer                        |
|          3 | isoform  | Transcript isoform ID within the gene.                 |

When some capture groups in a regular expression are named, and others are
unnamed, the unnamed capture groups are assumed to be the remaining capture
groups from the table above, in order by their position.

Any extra capture groups beyond those in the table above are ignored.

#### Examples

In the regular expression below, the default order of capture groups is
used. The first capture group (after the string `foo`, but before `bar`) is the
coverage. The second capture group (after the string `bar`, but before `baz`) is
the gene ID. The third and last capture group (after the string `baz`, but
before `_`) is the isoform ID.

```text
foo(.*)bar(.*)baz(.*)_.*
```

The regular expression below is the same, but capture groups are identified
explicitly.

```text
foo(?Pcoverage.*)bar(?Pgene.*)baz(?Pisoform.*)_.*
```

In the regular expression below, any text after the isoform ID is captured, but
RNA-clique will ignore it.

```text
foo(.*)bar(.*)baz(.*)_(.*)
```

In the following regular expression, the first capture group is the gene ID. The
second represents the coverage, and the third represents the isoform ID.

```text
tane(?Pgene.*)slafto(.*)dogge(.*)
```

### keep\_all

The last step in creating a [gene matches table](formats.md#gene-matches-tables)
is selecting the top gene pair for each sample 1 gene by bitscore. When
`keep_all` is False, this step produces a table such that every sample 1 gene in
the table is implicitly mapped to a single best match in sample 2. If for some
sample 1 gene there are multiple gene pairs with highest bitscore, ties are
broken arbitrarily by keeping only the row that comes first in the table.

When `keep_all` is True, RNA-clique allows more than one gene pair to be kept
for a sample 1 gene in the case of ties.

### path\_to\_sample

The `path_to_sample` setting should be a `dict` (YAML mapping) mapping [top
genes](formats.md#top-genes) files to the names of the sames they represent. By
default, RNA-clique assumes that each sample are named after the directory in
which its traanscripts FASTA file is found. This default behavior is reflected
in output YAML files produced by RNA-clique, but the mapping can be changed
manually to assign different names to samples.

#### Examples

The example below uses the set of 16 tall fescue samples analyzed in the
original RNA-clique paper. All top genes files here are mapped to their default
sample names.

```yaml
path_to_sample:
  f16_rna_clique_out/od1/SRR2321383_top.fasta: SRR2321383
  f16_rna_clique_out/od1/SRR2321384_top.fasta: SRR2321384
  f16_rna_clique_out/od1/SRR2321385_top.fasta: SRR2321385
  f16_rna_clique_out/od1/SRR2321386_top.fasta: SRR2321386
  f16_rna_clique_out/od1/SRR2321387_top.fasta: SRR2321387
  f16_rna_clique_out/od1/SRR2321388_top.fasta: SRR2321388
  f16_rna_clique_out/od1/SRR7990321_top.fasta: SRR7990321
  f16_rna_clique_out/od1/SRR7990322_top.fasta: SRR7990322
  f16_rna_clique_out/od1/SRR8003736_top.fasta: SRR8003736
  f16_rna_clique_out/od1/SRR8003737_top.fasta: SRR8003737
  f16_rna_clique_out/od1/SRR8003753_top.fasta: SRR8003753
  f16_rna_clique_out/od1/SRR8003754_top.fasta: SRR8003754
  f16_rna_clique_out/od1/SRR8003755_top.fasta: SRR8003755
  f16_rna_clique_out/od1/SRR8003756_top.fasta: SRR8003756
  f16_rna_clique_out/od1/SRR8003761_top.fasta: SRR8003761
  f16_rna_clique_out/od1/SRR8003762_top.fasta: SRR8003762
```

In the example below, the samples have been renamed to reflect their genotypes.

```yaml
path_to_sample:
  f16_rna_clique_out/od1/SRR2321383_top.fasta: SRR2321383
  f16_rna_clique_out/od1/SRR2321384_top.fasta: SRR2321384
  f16_rna_clique_out/od1/SRR2321385_top.fasta: SRR2321385
  f16_rna_clique_out/od1/SRR2321386_top.fasta: SRR2321386
  f16_rna_clique_out/od1/SRR2321387_top.fasta: SRR2321387
  f16_rna_clique_out/od1/SRR2321388_top.fasta: SRR2321388
  f16_rna_clique_out/od1/SRR7990321_top.fasta: SRR7990321
  f16_rna_clique_out/od1/SRR7990322_top.fasta: SRR7990322
  f16_rna_clique_out/od1/SRR8003736_top.fasta: SRR8003736
  f16_rna_clique_out/od1/SRR8003737_top.fasta: SRR8003737
  f16_rna_clique_out/od1/SRR8003753_top.fasta: SRR8003753
  f16_rna_clique_out/od1/SRR8003754_top.fasta: SRR8003754
  f16_rna_clique_out/od1/SRR8003755_top.fasta: SRR8003755
  f16_rna_clique_out/od1/SRR8003756_top.fasta: SRR8003756
  f16_rna_clique_out/od1/SRR8003761_top.fasta: SRR8003761
  f16_rna_clique_out/od1/SRR8003762_top.fasta: SRR8003762
```

### matrix

The `matrix` setting is a path to the [distance
matrix](formats.md#distance-matrix).

### subset\_of

When one analysis uses a subset of samples from another, reusing its [top
genes](formats.md#top-genes) and [gene matches
tables](formats.md#gene-matches-tables), the `subset_of` setting should be set
in the subset configuration file to a file path to the superset configuration
file. 

The [`make_subset.py`](usage.md#make-subsetpy) script uses `subset_of` to
determine which analysis to subset, but no other scripts in RNA-clique currently
read and use the `subset_of` parameter. Nevertheless, `subset_of` is useful as
metadata to keep track of the relationships between analyses.

## Troubleshooting configuration

If RNA-clique doesn't seem to be working as you expected, or if you would simply
like to verify that RNA-clique is processing your command-line arguments and/or
configuration file correctly before running a long analysis, it can be helpful
to see exactly what parameters RNA-clique will use to run a script. All scripts
in RNA-clique that support configuration files provide this feature via the
`--show-config` option, which allows the user to see the configuration or
arguments that RNA-clique will use.

By default, `--show-config` prints to standard output a YAML configuration file
containing the effective settings that will be used by the program. Arguments
can also be provided along with `--show-config` to show the original "raw"
parsed command-line arguments (`original_args`), the processed command-line
arguments (`args`), or the processed configuration (`config`, the
default). Multiple such arguments can be provided to show multiple sets of
settings. For example, `--show-config args config` will show both the processed
command-line arguments and the processed configuration.

The format used for showing the requested sets of settings can be controlled
with the `--show-config-format` option. By default, YAML (`yaml`) is preferred
for displaying the processed configuration, but the sets of settings will
instead be displayed as a `dict` of Python `reprs`s (`dict`) by default when
`args` or `original_args` have been requested. In addition to the `yaml` and
`dict` formats, RNA-clique supports displaying sets of settings in JSON
(`json`). If `yaml` or `json` have been specified and the requested sets of
settings can not be serialized automatically in those formats, RNA-clique will
produce an error message.
