# RNA-clique Python API

This document describes some options for using RNA-clique in your own Python
code. Most scripts included in RNA-clique (described in the [Command-line usage
guide](../usage.md)) have corresponding functions that can be called to perform the
same operations of the script from Python code. Just as different scripts offer
finer-grained control of the RNA-clique analysis, the corresponding functions
allow for customized analyses from within custom Python code.

## Running a full RNA-clique analysis

The `rna_clique` function in the 
{{file_link("`rna_clique.py` script", "rna_clique.py")}}
provides a way to perform a full RNA-clique analysis from Python code.

Unlike the `rna_clique.py` script, the `rna_clique` function must be provided
with explicit paths for the outputs. The paths to provide are summarized in the
table below.

| Formal parameter | Description                                                                                                           |
|------------------|-----------------------------------------------------------------------------------------------------------------------|
| `out_dir_1`      | Directory for the top $n$ genes of each transcriptome.                                                                |
| `out_dir_2`      | Directory for gene matches tables.                                                                                    |
| `cache_dir`      | Directory for the `simple_blast` [BLAST database caches](https://github.com/actapia/simple_blast/README.md#db-caches) |
| `output_graph`   | Output gene matches graph.                                                                                            |
| `output_matrix`  | Output distance matrix.                                                                                               |

The caller must also provide the number of top genes to use via the `top_genes`
parameter.

```python
from rna_clique import rna_clique
from pathlib import Path

out_dir = Path("rna_clique_out")
out_dir.mkdir(exist_ok=True)
# Get the SampleSimilarity object and a dict mapping paths to their sample
# names.
sim, path_to_sample = rna_clique(
    [
        Path("path/to/transcriptome1_dir"),
        Path("path/to/transcriptome2_dir"),
        Path("path/to/transcriptome3_dir"),
    ],
	out_dir_1=out_dir / "od1",
	out_dir_2=out_dir / "od2",
	cache_dir=out_dir / "db_cache",
    output_graph=output_dir / "graph.pkl",
    output_matrix=output_dir / "matrix.h5",
	top_genes=50000
)
print(sim.get_dissimilarity_df())
```

The `rna_clique` function returns two values. The first value is a
`SampleSimilarity` object that can be used to get the "dissimilarity" (distance)
matrix as a [Pandas
DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html)
using the `get_similarity_df` method as in the example above. 

`SampleSimilarity` objects can also provide other information about the
filtering process.

```python
from multiset_key_dict import MultisetKeyDict

# Get number of samples.
number_of_samples: int = sim.sample_count

# Get a list of samples used in the analysis.
samples: list[Path] = sim.samples

# Get a DataFrame containing the genes found in ideal components.
ideal_component_genes: pd.DataFrame = sim.valid

# Get the gene matches tables.
tables: MultisetKeyDict = sim.comparison_dfs
# Get specific gene matches table.
transcriptome1_vs_2: pd.DataFrame = tables[
    [
	    "path/to/trancsriptome1_dir",
		"path/to/transcriptome2_dir"
	]
]

# Get gene matches tables, filtered to include only genes in ideal components.
filtered_tables = MultisetKeyDict(sim.restricted_comparison_dfs())

# Get gene matches graph.
graph: nx.Graph = sim.graph
```

## Run individual RNA-clique steps

It is possible to use the RNA-clique Python API to run individual steps of the
RNA-clique method. As of version 0.2, RNA-clique supports running all steps
through its Python API.

### Select top genes by k-mer coverage

You can select top $n$ genes by $k$-mer coverage and save the results for a
single transcriptome using the `select_top_and_save` function of the
`select_top_genes_all.py` module. `select_top_genes_and_save` takes at least
four parameters and returns the name of the output file and the inferred sample
name.

| Formal parameter      | Description                                                          | Default          |
|-----------------------|----------------------------------------------------------------------|------------------|
| `out_dir`             | Output directory in which to save top $n$ genes.                     |                  |
| `transcripts`         | Name (basename) of FASTA file containing transcripts.                |                  |
| `x`                   | Directory containing the transcripts FASTA file.                     |                  |
| `top`                 | Number of top genes to select.                                       |                  |
| `parse_transcript_id` | Function for [parsing transcript IDs](#working-with-transcript-ids). | `default_parser` |

`default_parser` in the table above is `transcripts.default_parser`, a function
that uses the default regular expression designed for rnaSPAdes assemblies,
`^.*cov_([0-9]+(?:\.[0-9]+))_g([0-9]+)_i([0-9]+)`. For more information about
functions for parsing transcript IDs, see the ["Working with transcript
IDs"](#working-with-transcript-ids) section.

```python
from select_top_genes_all import select_top_and_save

# Select top genes for "path/to/transcriptome1/transcripts.fasta" and save to
# "top_genes_dir"
top_genes_file, sample_name = select_top_genes_and_save(
    "top_genes_dir",
    "transcripts.fasta",
    Path("path/to/transcriptome1"),
	50000,
)
```

If you don't wish to save the selected top genes, you can also use the
`TopGeneSelector `class from the `select_top_genes.py` module. `TopGeneSelector`
can work on genes that are not read from a file; to accomplish this, it requires
as the first argument to its constructor a nullary function that produces an
iterator of `Bio.SeqRecord` objects from which to select top genes. If you would
rather simply pass `TopGeneSelector` a FASTA file from which to select top
genes, use the `from_path` classmethod instead of the constructor. `from_path`
requires at least two arguments.

| Formal parameter      | Description                                                          | Default          |
|-----------------------|----------------------------------------------------------------------|------------------|
| `path`                | Path to FASTA file containing transcripts.                           |                  |
| `top`                 | Number of top genes to select.                                       |                  |
| `parse_transcript_id` | Function for [parsing transcript IDs](#working-with-transcript-ids). | `default_parser` |

```python
from select_top_genes import TopGeneSelector
selector = TopGeneSelector.from_path(
    Path("path/to/transcriptome1/transcripts.fasta"),
	50000
)
```

Like `select_top_genes_and_save`

To get the top genes, use the `get_top_gene_seqs` method.

```python
top_genes: list[Bio.SeqRecord] = list(selector.get_top_gene_seqs())
```

Alternatively, to just get the gene IDs of the top genes (which is faster), use
the `get_top_genes` method instead.

```python
top_gene_ids: list[int] = list(select.get_top_genes())
```

### Getting gene matches tables

You can obtain gene matches tables for all pairs of some set of samples via the
`find_all_pairs` function from the `find_all_pairs.py` module. `find_all_pairs`
requires at least four arguments and returns an iterator over the gene matches
tables, an iterator over the paths to the gene matches tables, and the total
number of gene matches tables computed. The arguments to `find_all_pairs` are
summarized in the table below.

| Formal parameter | Description                                                                                                           | Default value                     |
|------------------|-----------------------------------------------------------------------------------------------------------------------|-----------------------------------|
| `inputs`         | Paths to top $n$ genes for samples.                                                                                   |                                   |
| `output_dir`     | Directory in which to save gene matches tables.                                                                       |                                   |
| `cache_dir`      | Directory for the `simple_blast` [BLAST database caches](https://github.com/actapia/simple_blast/README.md#db-caches) |                                   |
| `path_to_sample` | Function mapping transcript FASTA file paths to sample names.                                                         |                                   |
| `hf_args`        | Arguments to pass to `HomologFinder`.                                                                                 | `[]`                              |
| `jobs`           | Number of parallel jobs to use.                                                                                       | `multiprocessing.cpu_count() - 1` |

`find_all_pairs` uses the `find_homologs_and_save` function, which finds the
gene matches table for just a single pair of samples. `find_all_pairs` needs at
least three arguments and returns the gene matches table that was computed and
saved. The arguments for `find_homologs_and_save` are summarized below.

| Formal parameter | Description                                   | Default value |
|------------------|-----------------------------------------------|---------------|
| `transcripts1`   | Path to top $n$ genes of first sample.        |               |
| `transcripts2`   | Path to top $n$ genes of second sample.       |               |
| `out_path`       | File in which to save gene matches table.     |               |
| `hf_args`        | Arguments to pass to `HomologFinder`.         | `None`        |
| `hf_kwargs`      | Keyword arguments to pass to `HomologFilter`. | `None`        |

In turn, `find_homologs_and_save` uses the `HomologFinder` class from the
`find_homologs.py` module. `HomologFinder` can find the gene matches table for a
pair of samples without saving it to a file. Many of `HomologFinder`'s
constructor's arguments are relevant to the BLAST search used for getting a gene
matches table. The constructor arguments are described below.

| Formal parameter      | Description                                                          | Default value |
|-----------------------|----------------------------------------------------------------------|---------------|
| `parse_transcript_id` | Function for [parsing transcript IDs](#working-with-transcript-ids). |               |
| `top_n`               | Number of top hits to select for each query gene (big $N$).          |               |
| `evalue`              | e-value cutoff to use for BLAST search.                              |               |
| `keep_all`            | Keep all matches in the case of ties.                                |               |
| `debug`               | Enable debug behavior.                                               | `False`       |

If any additional keyword arguments are given to the `HomologFinder`
constructor, they are passed to the constructor for the `simple_blast` BLAST
search used.

Once a `HomologFinder` has been constructed, its `get_match_table` method can be
used to find the gene matches table for two sets of top $n$ genes. The paths to
these top $n$ genes are the parameters to the `get_match_table` method.

```python
from transcripts import default_parser

# Get all gene matches tables with find_all_pairs.
from find_all_pairs import find_all_pairs

tables_iter, paths_iter, count = find_all_pairs(
    [
        Path("top/transcriptome1_top.fasta"),
        Path("top/transcriptome2_top.fasta"),
        Path("top/transcriptome2_top.fasta"),
	],
	Path("gene_matches_tables"),
	Path("db_cache"),
	# Could also use just path_to_sample from path_to_sample.py here.
	path_to_sample={
        Path("top/transcriptome1_top.fasta"): "transcriptome1",
        Path("top/transcriptome2_top.fasta"): "transcriptome2",
        Path("top/transcriptome2_top.fasta"): "transcriptome3", 
    }.__getitem__,
    hf_args=[
        default_parser,
        1,
        1e-99,
        True
    ]
)

# Get one gene matches table with find_hommologs_and_save.
from find_all_pairs import find_homologs_and_save

table = find_homologs_and_save(
    Path("top/transcriptome1_top.fasta"),
    Path("top/transcriptome2_top.fasta"),
    Path("gene_matches_tables/transcriptome1--transcriptome2.h5"),
    hf_args=[
        default_parser,
        1,
        1e-99,
        True
    ],
)

# Get one gene matches table with FindHomologs.
from find_homologs import HomologFinder

finder = HomologFinder(
    default_parser,
	1,
	1e-99,
	True
)
tble = finder.get_match_table(
    Path("top/transcriptome1_top.fasta"),
    Path("top/transcriptome2_top.fasta"),
)
```

### Getting a gene matches graph

A gene matches graph can be constructed from gene matches tables using the
`build_graph` function in the `build_graph.py` module. `build_graph` takes only
one argument&mdash;an iterable of Pandas DataFrames representing the gene
matches tables.

For information on loading gene matches tables that have already been saved to
disk, see the ["Loading gene matches tables"](#loading-gene-matches-tables)
section.

```python
from build_graph import build_graph
graph = build_graph(tables)
```

### Getting ideal components

You can obtain the ideal components from a gene matches graph using the
`get_ideal_components` function in the `filtered_distance.py` module. 

`get_ideal_components` requires two arguments and returns a list of NetworkX
`Graph` objects representing the ideal components. The arguments required by
`get_ideal_components` are summarized in the table below.

| Formal parameter | Description                                   |
|------------------|-----------------------------------------------|
| `g`              | The graph from which to get ideal components. |
| `samples`        | The number of samples in the analysis.        |


```python
from filtered_distance import get_ideal_components

ideal_components: list[nx.Graph] = get_ideal_components(
    graph,
	number_of_samples
)
```

### Filtering gene matches tables

Gene matches tables are most easily filtered to include only genes from ideal
components using the `SampleSimilarity` class from the `filtered_distance.py`
module. `SampleSimilarity`'s constructor arguments are summarized below.

| Formal parameter | Description                        |        |
|------------------|------------------------------------|--------|
| `graph`          | Gene matches graph                 |        |
| `comparison_dfs` | Iterable of gene matches tables.   |        |
| `sample_count`   | Number of samples in the analysis. | `None` |

When `sample_count` is not provided, `SampleSimilarity` counts the number of
samples from the gene matches graph.

Once a `SampleSimilarity` object has been constructed, the filtered gene matches
tables can be retrieved via the `restricted_comparison_dfs` method, which yields
pairs. The first element in each pair is a `frozenset` containing the samples
that the table is for. The second element is the table itself. The iterator
returned by `restricted_comparison_dfs` can be passed to `MultisetKeyDict`'s
constructor to get a mapping from pairs of samples to filtered gene matches
tables.

```python
from filtered_distance import SampleSimilarity
from multiset_key_dict import MultisetKeyDict

sim = SampleSimilarity(graph, tables)
filtered_tables = MultisetKeyDict(sim.restricted_comparison_dfs())

# Get an individual table.
table = filtered_tables[
    [
	    Path("od1/transcriptome1_top.fasta"),
		Path("od1/transcriptome2_top.fasta")
	]
]
```

### Getting a genetic distance

You can obtain *similarities* from filtered gene matches tables using the
`similarities_from_dfs` function from the `similarity_computer.py`
module. `similarities_from_dfs` accepts just one argument&mdash;an iterable of
pairs. The first element of each pair should be an unordered pair of samples
that the filtered table is for. The second element of each pair should be the
table itself.

```python
from similarity_computer import similarities_from_dfs

similarities = similarities_from_dfs(filtered_tables)
```

Each distance can then be computed as `1 - similarity`.

```python
from multiset_key_dict import MultisetKeyDict
distances = MultisetKeyDict((k, 1 - v) for (k, v) in similarities)
```

## Working with transcript IDs

RNA-clique expects to be able to read various metadata about transcripts in
[transcriptomes](../formats.md#transcriptomes) from the transcripts' FASTA
headers (also sometimes called "transcript IDs" in this guide). In order to
retrieve such metadata from a transcript FASTA header, RNA-clique must be able
to **parse** the transcript IDs.

RNA-clique was originally designed with rnaSPAdes assemblies in mind, and
although the defaults in both the command line and Python API reflect this
original design, RNA-clique is flexible enough to allow custom parsing for
non-SPAdes assemblies. At the command-line, RNA-clique provides this flexibility
by allowing the user to provide a [custom regular
expression](../config.md#transcript_id_regex) for parsing transcript IDs. In its
Python API, RNA-clique provides even more flexibility via functional
programming.

A transcript ID parser is a unary function that maps a transcript's FASTA ID to
a `TranscriptID` object. A `TranscriptID` object is a kind of `namedtuple` with
three fields, summarized in the table below.

| Position | Type    | Name     | Description                                            |
|----------|---------|----------|--------------------------------------------------------|
| 0        | `float` | coverage | $k$-mer coverage, expressed as a floating-point number |
| 1        | `int`   | gene     | Gene ID, a non-negative integer                        |
| 2        | `int`   | isoform  | Transcript isoform ID within the gene.                 |

Unlike an ordinary `namedtuple`, `TranscriptID` automatically casts the
parameters at its construction to the types shown in the table, ensuring that
the values of the fields always have the types given.

For example, a function that simply splits the FASTA ID by colons (`:`) and
takes the last three elements as the transcript ID is shown below.

```python
from transcripts import TranscriptID

def custom_transcript_id_parser(x):
    return TranscriptID(*x.split(":")[-3:])
```

You can also create a parser from a regular expression in the same fashion as
the command-line interface using the `TranscriptID.parser_from_re`
classmethod. For example,

```python
import re

my_parser = TranscriptID.parser_from_re(re.compile(r"([^:]*):([^:]*):([^:]*)$"))
parsed = my_parser("12:13:14")
print(parsed) # TranscriptID(coverage=12.0, gene=13, isoform=14)
```

## Loading gene matches tables

An individual [gene matches table](../formats.md#gene-matches-tables) file can
be loaded using the `read_table` function from
`gene_matches_tables.py`. `read_table` requires at least one argument, the path
to the gene matches table to load, and it returns a Pandas DataFrame
representing the loaded gene matches table.

```python
from gene_matches_tables import read_table

table = read_table(Path("od2/transcriptome1--transcriptome2.h5"))
```

To get a list of gene matches tables in a directory, use `get_table_files` from
the same module.

```python
from gene_matches_tables import get_table_files

table_files = get_table_files(Path("od2"))
tables = [read_table(t) for t in table_files]
```

## Exporting orthologs (ideal components)

You can export orthologs from Python code using the `OrthologExporter` class
found in the `export_orthologs.py` module. The `OrthologExporter` class requires
at least two arguments for construction. The parameters accepted by the
constructor are summarized in the table below.

| Formal parameter      | Description                                                                     | Default value |
|-----------------------|---------------------------------------------------------------------------------|---------------|
| `sim`                 | `SampleSimilarity` for an analysis.                                             |               |
| `parse_transcript_id` | Function for [parsing transcript IDs](#working-with-transcript-ids).            |               |
| `non_contributing`    | Whether to filter components where there are no differences in aligned regions. | `True`        |
| `consistent_strands`  | Whether to reorient transcripts to be consistent per ideal component.           | `True`        |
| `allow_inconsistent`  | Whether to continue when basic transcript reorientation method fails.           | `True`        |
| `jobs`                | Number of parallel jobs to use.                                                 | `1`           |
| `debug`               | Whether to enable debug behavior.                                               | `False`       |

After construction, calling an `OrthologExporter` object's `by_component` or
`by_sample` method with export the orthologs, organizing them into files by
component or by sample, respectively. The `by_sample` and `by_component` methods
also accept a few different options. These options are summarized in the table
below.

| Formal parameter | Description                                                                           | Default value |
|------------------|---------------------------------------------------------------------------------------|---------------|
| `out_dir`        | Output directory for exported orthologs.                                              |               |
| `rename`         | Function to give a new name to a transcript from its sample ID, gene, and isoform ID. | `None`        |
| `order`          | Whether to put the new name before or after the original name.                        | `after`       |
| `set_rlimit`     | (`by_component` only) Try to Set `rlimit` for writing many files at once.             | `True`        |
| `make_all`       | Create combined file containing all exported orthologs.                               | `True`        |

Importantly `make_all` **must be set to `True`** if you wish to later [search
the ortholog sequences](#searching-ideal-components-exported-orthologs) because
the `search` function requires a path to the combined file with all the
exported orthologs. 

`by_sample` and `by_component` both return dictionaries mapping each groups by
which the orthologs were organized to the paths to the file in which the
orthologs in that group were saved. For `by_sample`, the `dict` maps each sample
name to the file in which orthologs from that sample are stored. For
`by_component`, the `dict` maps each ideal component ID (expressed as an
integer) to the file in which orthologs from that ideal component are stored.

```python
from export_orthologs import OrthologExporter
from transcripts import default_parser

exporter = OrthologExporter(sim, default_parser)
component_to_path = exporter.by_component(Path("exported_orthologs"))
```

## Searching ideal components (exported orthologs)

Exported orthologs can be searched using the `search` function from the
`search_ideal_components.py` module. `search` requires at least five arguments
but can accept many optional arguments. The possible parameters are described in
the table below.

| Formal parameter      | Description                                                                                                            | Default value    |
|-----------------------|------------------------------------------------------------------------------------------------------------------------|------------------|
| `sim`                 | `SampleSimilarity` for analysis in which search is performed.                                                          |                  |
| `exported`            | Path to *combined* FASTA file containing all exported orthologs.                                                       |                  |
| `db_cache_loc`        | Directory for the `simple_blast` [BLAST database caches](https://github.com/actapia/simple_blast/README.md#db-caches). |                  |
| `out_dir`             | Output directory in which to store search results.                                                                     |                  |
| `query`               | Path to FASTA file containing query sequences.                                                                         |                  |
| `parse_transcript_id` | Function for [parsing transcript IDs](#working-with-transcript-ids).                                                   | `default_parser` |
| `path_to_sample`      | Function retrieving sample names from top gene paths.                                                                  | `path_to_sample` |
| `extended_evalue`     | e-value cutoff for extended search. `None` disables extended search.                                                   | `None`           |
| `export_components`   | Whether to export components in which matchse are found.                                                               | `True`           |
| `merge_sams`          | Whether to merge extended searach results into one SAM file.                                                           | `False`          |
| `strand_graph`        | Orientation graph (strand graph) to use. Will be created if not provided.                                              | `False`          |
| `node_to_ccc`         | Mapping from gene ID, isoform ID pairs to components in the meta-strand graph. Will be created if not provided.        | `None`           |
| `evalue`              | e-value cutoff to use for initial searches.                                                                            | `None`           |
| `jobs`                | Number of parallel jobs to use.                                                                                        | `1e-50`          |
| `debug`               | Enable debug behavior.                                                                                                 | `False`          |

If you have an `OrthologExporter` object, you can speed `search` up by providing
it with the `OrthologExporter`'s `strand_graph` and
`node_to_component_component` attributes for its `strand_graph` and
`node_to_ccc` parameters, respectively,

`search` returns a `SearchResult` object, which is a `namedtuple` containing
attributes `hits`, `seqs`, and `components`. `hits` is the number of BLAST
HSPs to the ideal components. `seqs` is the number of transcripts producing
HSPs. `components` is the number of components producing HSPs.

```python
from search_ideal_components import search, SearchResult

result: SearchResult = search(
    sim,
	Path("exported_orthologs/all_ideal.fasta"),
	Path("db_cache"),
	Path("search_results),
	Path("queries.fasta"),
)
# Print the number of total HSPs, transcripts and components producing HSPs.
print(*result)

# Get a faster result by providing strand_graph and node_to_cc from the
# OrthologExporter.
fast_result: SearchResult = search(
    sim,
	Path("exported_orthologs/all_ideal.fasta"),
	Path("db_cache"),
	Path("search_results),
	Path("queries.fasta"),
	strand_graph=exporter.strand_graph,
	node_to_ccc=exporter.node_to_component_component,
)
```

## Working with configuration objects

RNA-clique supports configuration via [YAML configuration files](../config.md)
as well as via a command-line interface. Additionally, some RNA-clique programs
will automatically produce a configuration file for reproducing an
analysis. Althoug using RNA-clique configuration files in not necessary to use
the Python API, it can sometimes be convenient to be able to read or create
configuration files from your own Python code.

To load an RNA-clique configuration file, use the `RNACliqueConfig.yaml_load`
classmethod from the `config.py` module. `yaml_load` takes just one
argument&mdash;the path to the configuration file&mdash;and loads and returns
the configuration at that path.

```python
import config as config_module

config = config_module.RNACliqueConfig.yaml_load(Path("path/to/config.yaml"))
```

The attributes of configuration objects are documented in the [Configuration
guide](../config.md). To save an `RNACliqueConfig` object, use the `yaml_save`
method. `yaml_save` takes one argument, the output path.

```python
config.yaml_save(Path("path/to/config.yaml"))
```

## Making subsets

To create a subset of an existing analysis from Python code, you can use the
`SubsetAnalysisCreator` class from the `make_subset.py`
module. `SubsetAnalysisCreator`'s constructor requires three arguments, which
are summarized in the table below.

| Formal parameter | Description                                                                                            |
|------------------|--------------------------------------------------------------------------------------------------------|
| `matches`        | Predicate indicating whether a path should be included in the subset.                                  |
| `super_config`   | [`RNACliqueConfig`](#working-with-configuration-objects) for the parent (original, superset) analysis. |
| `config`         | Partial [`RNACliqueConfig`](#working-with-configuration-objects) for the subset analysis.              |

Each of these parameters requires some further explanation.

First, the `matches` parameter is a predicate&mdash;i.e., a function that
returns a `bool`. It should take the `Path` to a [top
genes](../formats.md#top-genes) file and return `True` if the `Path` should be
included in the subset, or `False` otherwise. The accepted `matches` parameter
is a function to provide maximum flexibility, but some other functions
(discussed below) are also provided to make creating such functions easier.

Second, the `super_config` parameter is the `RNACliqueConfig` for the original
analysis to be subsetted. The `super_config` is used to get the location of the
original gene matches tables and the paths to the top genes files.

Third, the `config` parameter is the `RNACliqueConfig` for the subset
analysis. Since the subset has not yet been created, `config` is not expected to
have all of its attributes assigned. The only required (non-`None`) attributes
are `table_dir` and `graph`&mdash;these paths are where the
`SubsetAnalysisCreator` will put the gene matches table links and new gene
matches graph for the subset analysis. `SubsetAnalysisCreator` does not modify
either `RNACliqueConfig` passed to it in the constructor but does create and
modify a copy of `config`.

After constructing a `SubsetAnalysisCreator`, calling the object's `make` method
creates the analysis using the parameters provided to the constructor. After
`make` has been called, the `config` attribute of the `SubsetAnalysisCreator`
instance will be a complete `RNACliqueConfig` for the subset analysis.


```python
from make_subset import SubsetAnalysisCreator

child_config = RNACliqueConfig()
child_config.tables_dir = Path("subset/tables_dir")
child_confid.graph = Path("subset/graph.pkl")
subsetter = SubsetAnalysisCreator(
    # Exclude samples that start with the letter "B".
    lambda x: not x.name.startswith("B"),
	parent_config,
	child_config
)
```

Since creating a child analysis `RNACliqueConfig` just to provide it to
`SubsetAnalysisCreator` is somewhat cumbersome, `SubsetAnalysisCreator` also
provides a classmethod, `from_paths` that allows one to simply provide the
subset analysis' `tables_dir` and `graph` directly to
`SubsetAnalysisCreator`. The parameters accepted by that classmethod are
described in the table below.

| Formal parameter | Description                                                                                            |
|------------------|--------------------------------------------------------------------------------------------------------|
| `matches`        | Predicate indicating whether a path should be included in the subset.                                  |
| `super_config`   | [`RNACliqueConfig`](#working-with-configuration-objects) for the parent (original, superset) analysis. |
| `tables_dir`     | Directory in which to store child's links to parent's gene matches tables.                             |
| `graph`          | Path at which to create gene matches graph for child analysis.                                         |


```python
from make_subset import SubsetAnalysisCreator

subsetter = SubsetAnalysisCreator.from_paths(
    # Exclude samples that start with the letter "B".
    lambda x: not x.name.startswith("B"),
	parent_config,
    Path("subset/tables_dir"),
	Path("subset/graph.pkl"),	
)
```

## Visualization

Although visualization is not one of the primary features of RNA-clique,
RNA-clique nevertheless provides some code designed for visualization. The API
for those function is described [in a separate document](visualization.md).
