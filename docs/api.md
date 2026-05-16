# RNA-clique Python API

This document describes some options for using RNA-clique in your own Python
code. Most scripts included in RNA-clique (described in the [Command-line usage
guide](usage.md)) have corresponding functions that can be called to perform the
same operations of the script from Python code. Just as different scripts offer
finer-grained control of the RNA-clique analysis, the corresponding functions
allow for customized analyses from within custom Python code.

## Running a full RNA-clique analysis

The `rna_clique` function in the [`rna_clique.py` script](/rna_clique.py)
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
[transcriptomes](formats.md#transcriptomes) from the transcripts' FASTA headers
(also sometimes called "transcript IDs" in this guide). In order to retrieve
such metadata from a transcript FASTA header, RNA-clique must be able to
**parse** the transcript IDs.

RNA-clique was originally designed with rnaSPAdes assemblies in mind, and
although the defaults in both the command line and Python API reflect this
original design, RNA-clique is flexible enough to allow custom parsing for
non-SPAdes assemblies. At the command-line, RNA-clique provides this flexibility
by allowing the user to provide a [custom regular
expression](config.md#transcript_id_regex) for parsing transcript IDs. In its
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

An individual [gene matches table](formats.md#gene-matches-tables) file can be
loaded using the `read_table` function from
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

RNA-clique supports configuration via [YAML configuration files](config.md) as
well as via a command-line interface. Additionally, some RNA-clique programs
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
guide](config.md). To save an `RNACliqueConfig` object, use the `yaml_save`
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
genes](formats.md#top-genes) file and return `True` if the `Path` should be
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
	
## Visualizations

RNA-clique does not primarily aim to offer tools for data visualization; other
software exists for that and probably does it better than RNA-clique
could. Nevertheless, some data visualization functions that the author uses for
his own analyses are included with RNA-clique and can be used via the Python
API.

Examples for using these visualization functions are also provided as part of
the end-to-end ["From RNA-seq reads to a phylogenetic tree with
RNA-clique"](tutorials/reads2tree/README.md) tutorial. Sample code from these
tutorials is linked in the subsections covering individual visualizations below.

### Sample metadata

The visualization functions depend on sample metadata to properly organize,
label, and visually encode information about each sample in the created
visualizations. In this context, sample metadata should be a Pandas DataFrame in
which there is one-to-one correspondence between rows and samples. Each row
should give some information about the sample. For example, sample metadata
might include information about which population a sample comes from, or what
the sample's chemotype is.

It must be possible to map each sample to its corresponding metadata by the ID
of the sample in the distance matrix. To make this possible, one column should
contain these IDs. Which column to use can usually be specified to visualization
function via the `sample_name_column` keyword argument. Note that, by default,
`SampleSimilarity` gives distance matrices with rows and columns labeled with
the paths to the top genes for the sample rather than the sample name itself. If
your metadata uses the names of the samples, it may be necessary to rename the
rows and columns of the distance matrix to use the actual names of the samples.

### Heatmaps

A distance matrix can be visualized in a straightforward fashion via a
heatmap. In a heatmap, the elements of the distance matrix are represented by
rectangles in a grid. Labels for the rows and columns of the matrix are shown
above and to the left of the rows and columns of squares in the heatmap. Each
rectangle representing a matrix element is colored according to a colormap,
which is also shown with a scale next to the heatmap.

RNA-clique offers a `draw_heatmap` function, which can be found in the
`heatmap.py` module. The arguments accepted by `draw_heatmap` are summarized in
the table below.

| Formal parameter     | Description                                                          | Default value               |
|----------------------|----------------------------------------------------------------------|-----------------------------|
| `mat`                | The distance (or similarity) matrix to visualize as a heatmap.       |                             |
| `sample_metadata`    | Sample metadata for sorting and grouping samples in the heatmap.     | `None`                      |
| `sample_name_column` | Column of the sample metadata that names each sample.                | `"name"`                    |
| `order_by`           | Column(s) by which to order samples along the axes of the heatmap.   | `None`                      |
| `square`             | Whether to draw heatmap cells as squares.                            | `True`                      |
| `group_by`           | Column(s) by which to group samples along the axes of the heatmap.   | `None`                      |
| `draw_group_labels`  | Whether to draw labels for sample groupings.                         | `False`                     |
| `make_group_label`   | Function to get label for group from group value.                    | `default_group_label_maker` |
| `digit_annot`        | Digits of actual value to show within cells, if any.                 | `None`                      |
| `label_padding_x`    | Horizontal padding to add to labels on vertical axis.                | `0.0275`                    |
| `label_padding_y`    | Vertical padding to add to labels on horizontal axis.                | `0.0275`                    |
| `sort_key`           | Function to transform sort column(s).                                | `None`                      |
| `label_kwargs`       | Additional keyword arguments to pass to `plt.text` for group labels. | `None`                      |
| `x_label_kwargs`     | Additional keyword arguments for labels on vertical axis.            | `None`                      |
| `y_label_kwargs`     | Additional keyword arguments for labels on horizontal axis.          | `None`                      |
| `draw_debug_points`  | Whether to draw points useful for debugging drawing code.            | `False`                     |

Any additional arguments or keyword arguments provided to `draw_heatmap` will be
passed to Seaborn's
[`heatmap`](https://seaborn.pydata.org/generated/seaborn.heatmap.html) function,
which is what `draw_heatmap` uses as a base for its heatmap plots.

RNA-clique's `draw_heatmap` function assumes that the distance matrix to be
drawn is symmetric, so it only draws the upper triangle of the matrix,
regardless of whether the matrix it is given is actually symmetric.


Samples can be both ordered and grouped by values from the `sample_metadata`.
When samples are grouped by some columns, those samples that have the same value
for those columns will be adjacent in the rows and columns of the heatmap. When
samples are ordered by some columns, the samples are sorted
according to their values for those columns along the rows and columns of the
heatmap. If samples are ordered by some columns, then they are also grouped by
those columns, but not vice versa. Hence, the `group_by` columns must always be
a prefix of the `order_by` columns. If `order_by` is specified, but `group_by`
is not, then `group_by` will be set automatically to `order_by`.

Vertical and horizontal lines separate groups on the axes of the heatmap. When
`draw_group_label` is `True`, labels are also drawn to identify the groups. What
the labels say can be controlled via the `make_group_label`
function. `make_group_label` should take the value(s) for a particular group and
return a string to be used as the label for that group. By default, the group
label is simply the values (converted to strings) and joined with commas and
spaces.

As a supplement to the colormapping, it is sometimes helpful to display digits
of the actual value of a matrix element on top of the heatmap cell for that
element. Since distance values do not typically have short decimal
representations, only a few digits can usually be shown per heatmap cell. The
`digit_annot` parameter, when specified (not `None`), shows no more than the
specified number of digits of each value. The value shown in a cell represents
the value of the element in the smallest `1/(10**n)` unit such that the largest
value to show needs no more than `digit_annot` digits. Values shown are rounded
to the nearest unit. For example, if the distances to show are `0.123`, `0.012`,
and `0.015`, displaying with two digits would display values in hundredths,
giving `12`, `1`, and `2`, respectively.

In some cases, it may be desirable to sort the samples by values that don't
appear explicitly in the `sample_metadata`. The `draw_heatmap` function supports
such sorting via the `sort_key` keyword argument. `sort_key` must be a function
that maps the `order_by` columns in `sample_metadata` to other columns to be
used for sorting.

```python
from heatmap import draw_heatmap
from path_to_sample import path_to_sample

sample_metadata: pd.DataFrame = pd.read_csv("my_metadata.csv")
distances: pd.DataFrame = sim.get_dissimilarity_df().rename(
    columns=path_to_sample
).rename(
    rows=path_to_sample
)

# Draw a heatmap with samples grouped and ordered by population and chemotype.
# Display two digit representations of the distances on the cells.
draw_heatmap(
    distances,
	sample_metadata,
	order_by=["population", "chemotype"],
	draw_group_labels=True,
	digit_annot=2,
)
```

A real example of using `draw_heatmap` can also be found at
[`docs/tutorials/reads2tree/make_heatmap.py`](docs/tutorials/reads2tree/make_heatmap.py).

### PCoA plots

A Principal Coordinate Analysis (PCoA) plot represents samples as points in a
scatter plot. Point coordinates are chosen for the samples such that distances
are preserved as well as possible given the number of dimensions allowed.

RNA-clique provides the `draw_pcoa` function for drawing PCoA plots in two or
three dimensions. `draw_pcoa` is in the `pcoa.py` module and returns a
[`scikit-bio`](https://scikit.bio/) `skbio.stats.ordination.OrdinationResults`
object representing the results of the PCoA analysis. The parameters accepted by
`draw_pcoa` are summarized in the table below.

| Formal parameter        | Description                                                                                          | Default value               |
|-------------------------|------------------------------------------------------------------------------------------------------|-----------------------------|
| `dis_df`                | Distance matrix to visualize, as a Pandas DataFrame.                                                 |                             |
| `sample_metadata`       | Sample metadata to use for encoding sample information into the point markers and ellipsoids.        |                             |
| `group_by`              | Columns from sample metadata to encode in point markers and ellipsoids.                              |                             |
| `sample_name_column`    | Column containing sample names in `sample_metadata`.                                                 | `"name"`                    |
| `make_group_label`      | Functions mapping group values to text labels to display in legend.                                  | `default_group_label_maker` |
| `labelers`              | Functions to transform group values before making group labels.                                      | `None`                      |
| `colors`                | Colors to use for encoding values of the `group_by` column.                                          | `None`                      |
| `legend`                | Whether to create a legend.                                                                          | `True`                      |
| `index_to_kwargs`       | Function mapping index among `group_by` values to keyword arguments to pass to `scatter`.            | `empty_dict`                |
| `group_to_kwargs`       | Function mapping `group_by` values to keyword arguments to pass to `scatter`.                        | `empty_dict`                |
| `index_group_to_kwargs` | Function mapping both index and values for `group_by` to keyword arguments to pass to `scatter`.     | `empty_dict`                |
| `order_by`              | Columns on which to order samples (influences order that groups appear in legend).                   | `None`                      |
| `sort_key`              | Functions to transform `order_by` columns for purposes of sorting.                                   | `None`                      |
| `ellipsoids`            | Whether to draw ellipsoids on the plot.                                                              | `False`                     |
| `make_ellipsoid`        | Function to create ellipsoid to draw from data for a group.                                          | `conf_ellipsoid(0.95)`      |
| `ellipsoid_kwargs`      | Extra keyword arguments to pass to `draw_ellipsoid`.                                                 | `None`                      |
| `contribution`          | Display relative contributions of principal component axes on axis labels.                           | `True`                      |
| `annotate`              | Label individual points with their sample names.                                                     | `False`                     |
| `adjust`                | Whether to space out labels automatically using [`adjustText`](https://github.com/Phlya/adjustText). | `True`                      |
| `legend_factors`        | Try to separate independent encodings in the legend.                                                 | `False`                     |
| `default_legend_marker` | Default kwargs to pass to `mpl.lines.Line2D` for creating legend entries.                            | `default_marker_style`      |
| `dropna`                | Drop rows where the `group_by` columns have NA values.                                               | `True`                      |
| `dimensions`            | Number of dimensions for the PCoA plot.                                                              | `2`                         |
| `ax`                    | `matplotlib.axes.Axes` on which to draw the plot.                                                    | `None`                      |
| `axis_label_kwargs`     | Keyword arguments for creating axis labels.                                                          | `None`                      |

Although a few of these options are self-explanatory, many of them require
further explanation.

#### group\_by

When making a PCoA plot, it is common to want to encode metadata about samples
using various visual variables such as point marker color, shape, opacity,
etc. For example, one might want to encode population with marker shape and
chemotype with marker color. `draw_pcoa` enables encoding metadata in visual
variables by allowing the caller of the function to specify metadata columns on
which to group samples via the `group_by` parameter. All samples that have the
same value for all of the columns in `group_by` form a group. Other parameters
to `draw_pcoa` can then be used to specify how groups should be drawn
differently.

#### make\_group\_label and labelers

The actual values for the `group_by` columns are not necessarily what you want
displayed in the legend entries corresponding to those values. To allow the
caller to specify arbitrary labels for each `group_by` columns value,
`draw_pcoa` allows the caller to specify `make_group_label` and `labelers`.

The `labelers` are functions such that `labelers[i]` is used to create the label
for the value of column `i` of `group_by`. Each element of `labelers` should
accept a value of column `i` of `group_by` and should return a string, the label
to use for that value. If a labeler doesn't exist for column `i`&mdash;i.e., if
a labeler is `None`, or if `i >= len(labelers)`&mdash;then `draw_pcoa` treats
that column as if its labeler were the identity function. 

Each `labeler` element operates on a value from a single column and returns
a label for that value in that column. The `make_group_label` function operates
on an iterable of labels for all of the columns in `group_by` and returns one
final label to use for the group with those individual column
labels. By default, the value for `make_group_label` is the
`default_group_label_maker`, which returns the individual labels, joined by
commas and spaces.

#### colors

The simplest way of specifying how values of the `group_by` columns should be
encoded via visual variables is by using the `colors` parameter. As the name
suggests, `colors` only allows assigning colors to groups, but it provides a
convenient and flexible way of doing so. `colors` can be a
`matplotlib.colors.ListedColormap`, a sequence of RGB tuples, or a mapping from
tuples of `group_by` column values to RGB tuples. In the first two cases,
`group_by` values are assigned colors sequentially. The first value for the
`group_by` columns gets the first color. The second value gets the second color,
and so on. In the third case, `group_by` column values are assigned explicitly
to specific colors via the mapping.

#### index\_to\_kwargs, group\_to\_kwargs, and index\_group\_to\_kwargs

A more flexible (but more cumbersome) way of assigning visual encodings to
`group_by` values is provided by the `index_to_kwargs`, `group_to_kwargs`, and
`index_group_to_kwargs` parameters. Each of these parameters is a function that
takes some parameter(s) identifying the `group_by` column values and returns
keyword arguments to pass to `matplotlib`'s `scatter` function for the samples
in that group. `index_to_kwargs`, `group_to_kwargs`, and `index_group_to_kwargs`
differ in what parameters they should accept.

`index_to_kwargs` should accept the index of the `group_by` column value as its
only parameter. `index_to_kwargs` allows for sequential mappings like those that
are possible by providing the `colors` parameter with a `ListedColormap`. For
example, `index_to_kwargs` might map the first value (index `0`) to 
`{"color": "blue", "marker": "o"}` and map the second value (index `1`) to 
`{"color": "red", "marker": "s"}` to make the first value blue circles and the
second value red squares.

`group_to_kwargs` should accept the *value* of the `group_by` column itself as
its only parameter. `group_to_kwargs` allows for explicit mappings like those
that are possible by providing the `colors` parameter with a `Mapping`. For
example, `group_to_kwargs` might map `("location1", "chemotype2")` to 
`{"color": "blue", "marker": "^"}`, map `("location2", "chemotype1")` to
`{"color": "red", "marker": "s"}`, and map anything else to 
`{"color": "gray", "marker": "o"}`.

Finally, `index_group_to_kwargs` should accept both the index and the value of
the `group_by` column value as its two parameters, in that
order. `index_group_to_kwargs` is useful when some visual encodings need to be
specified explicitly, but others can be specified by a simple sequences. For
example, you could map `("location1",)` and `("location2",)` to 
`{"color": "red"}` and `{"color": "blue"}`, respectively, and have all other
values map to colors corresponding to their index in a `ListedColormap`.

#### order\_by and sort\_key

When `legend_factors` in not specified, entries in the legend produced by
`draw_pcoa` appear in the order in which they were drawn on the scatter plot. By
default, this order is the sorted order of the `group_by` values, but the order
can be changed by providing the `order_by` parameter to `draw_pcoa`. When
`order_by` is provided, all samples are sorted by `order_by` before drawing,
causing the entries to appear in the order in which their `group_by` values
appear in the sorted DataFrame.  In addition to the columns from the sample
metadata, `order_by` can be used with the PCoA coordinates of samples. This
feature can be used, for example, to ensure that groups that appear in the
legend in the same vertical order they appear in the plot.

In some cases, it is useful to sort the dataframe by values that do not appear
in the columns of the sample metadata but can be computed from one or more
columns. In that case, the `sort_key` parameter can be provided. `sort_key`
should be a function taking the `order_by` columns and returning a new Pandas
Series object that will be used to sort the samples.

#### ellipsoids and make\_ellipsoid

In some cases, it is desirable to draw an ellipsoid for each `group_by`
group. This behavior can be enabled by providing `True` to the `ellipsoids`
parameter. By default, when `ellipsoids=True`, `draw_pcoa` will draw a 95%
confidence ellipsoid for the population mean of each group of samples, assuming
the population distribution is multivariate normal. Different ellipsoids can be
drawn by providing a value to the `make_ellipsoid` parameter. `make_ellipsoid`
should be a function that accepts a data matrix in which rows represent
individual samples (observations), and columns represent the principal component
dimensions. The function should return an `Ellipsoid` from the
`confidence_ellipsoid.py` module, which is defined by its center and vectors
representing its axes (one for each dimension). The `confidence_ellipsoid`
module also contains some functions for creating such `make_ellipsoid`
functions, including `get_multivariate_normal_density_ellipsoid`, which returns
an ellipsoid containing the specific probability density, and
`get_multivariate_normal_confidence_ellipsoid`, which returns a confidence
ellipsoid at the given level for the population mean.

#### contribution

When `contribution` is specified, the relative contribution of each principal
component axis is shown on its axis label. The relative contribution of an axis
is computed as the eigenvalue of that axis divided by the sum of eigenvalues of
all possible axes (including ones greater than the number of dimensions).

#### annotate and adjust

When the `annotate` parameter is specified, individual samples will be labeled
with their names from the `sample_name_column`. Initially, each label is placed
exactly on the point to which it refers, but since this can cause labels to
overlap, the text positions might need to be adjusted. Such adjust takes place
automatically via [`adjustText`](https://github.com/Phlya/adjustText) when the
`adjust` parameter is `True` (the default); setting `adjust=False` keeps the
original positions.

#### legend\_factors and default\_legend\_marker

In some cases, encodings of column might be independent in the sense that the
value of one `group_by` column completely determines one variable of the visual
encoding, and vice versa, and another column completely determines another
variable, and vice versa. For example, consider the `group_by` values have the
mapping shown in the table below,

| `group_by` value | Marker keyword arguments            |
|------------------|-------------------------------------|
| `('a', '1')`     | `{"color": "red", "shape": "o"}`    |
| `('b', '2')`     | `{"color": "blue", "shape": "s"}`   |
| `('c', '1')`     | `{"color": "yellow", "shape": "o"}` |
| `('a', '2')`     | `{"color": "red", "shape": "s")`    |
| `('b', '1')`     | `{"color": "blue", "shape": "o"}`   |
| `{'c', '2'}`     | `{"color": "yellow", "shape": "s"}` |

In the example above, the value for the first `group_by` column determines the
color. When the value is `'a'`, the color is red. When the value is `'b'`, the
color is blue, and when the value is `'c'`, the color is yellow. The value of
the second variable has no impact on the color. Likewise, when the value of the
second variable is `'1'`, the shape is a circle, and when the value is `'2'`,
the shape is a square. The value of the first variable has no impact on shape.

By default, `draw_pcoa` would make a legend entry for each of the six
combinations of color and shape seen among the samples, but, in this case,
because we have variables that are encoded independently, we express the legend
more succinctly with just five entries. The five entries would should that red
corresponds to `'a'`, blue to `'b'`, and yellow to `'c'`, and that circle
corresponds to `'1`' and square to `'2'`. 

`draw_pcoa` can try to detect such independent variables when the
`legend_factors` parameter is set to a list of tuples or strings or to `True`
(equivalent to `legend_factors` being set to `group_by`). `draw_pcoa` will then
display each visual variable's mapping separately. This is accomplished by first
computing what keyword arguments are common to all samples that share values for
each of the provided tuples of columns in `legend_factors`. For the example
above, if the columns were named `alpha` and `beta`, and `legend_factors` was
set to `[('alpha',), ('beta',)]`, then `draw_pcoa` would recognize that all
cases where `alpha = 'a'` are cases where `color` is mapped to red, 
`alpha = 'b'`, blue, and `alpha = 'c'`, yellow. Likewise, it would recognize
that all cases where `beta = '1'` are cases where `shape` is mapped to `o`, and
all cases where `beta = '2'` are cases where `shape` is mapped to `s`.

The legend is then drawn one set of columns in `legend_factors` at a time; each
value for each set of columns is shown next to a marker using the shared set of
keyword arguments common to all samples with that value for that set of
columns. For the example above, `'a'` would be shown next to a red marker, `'b'`
next to blue, and `'c'` next to yellow. Then, `'1`' would be shown next to a
circle, and `'2'` would be shown next to a square.

Since a marker must have some value for the unspecified visual variables,
default values for those variables can be specified via
`default_legend_marker`. For example, the red marker shown in the example above
could not simply be a red marker&mdash;it would also have to have a shape,
opacity, etc. `default_legend_marker` allows you to specify those default
variables.

Since `draw_pcoa`'s behavior when `legend_factors=True` is simply based on
analyzing the provided mappings from `group_by` values to `scatter` keyword
arguments, it may sometimes give unexpected results, splitting variables that
are only coincidentally independent, or failing to split variables when certain
combinations of values are always observed together.

```python
from pcoa import draw_pcoa

# Create a PCoA plot ...
res: skbio.stats.ordination.OrdinationResults = draw_pcoa(
   sim.get_dissimilarity_df(),
   sample_metadata,
   # Group samples by location and chemotype.
   ["location", "chemotype"],
   "sample_name",
   # Samples in place1 with chemotype lolC- should be red. Everything else just
   # gets its color from the tab10 colormap.
   index_group_to_kwargs=lambda i, x: {
           ("place1", "lolC-"): {"color": "red"}
   }.get(
		x,
		{"color": matplotlib.colormaps.get("tab10").colors[i]}
   ),
   # Draw 95% confidence ellipsoids for population mean.
   ellipsoids=True
)
```

A real example of using `draw_pcoa` can also be found at
[`docs/tutorials/reads2tree/make_pcoa.py`](docs/tutorials/reads2tree/make_pcoa.py).

### Phylograms

RNA-clique provides functions for drawing phylograms in the `phylo_utils.py`
module. Unlike the functions for heatmaps and PCoA plots, the functions for
phylograms expect to be given a BioPython Tree instead of a distance
matrix. Such trees can be constructed using BioPython or some other Python
libraries. 

The main function for drawing phylograms is the `draw_tree` function in the
`phylo_utils.py` module. The parameters that the function accepts are described
in the table below.

| Formal parameter     | Description                                          | Default value |
|----------------------|------------------------------------------------------|---------------|
| `tree`               | The phylogenetic tree for which to draw a phylogram. |               |
| `blank_nonterminals` | Whether to set non-terminal labels to empty string.  | `True`        |
| `clades`             | Mapping from values to clades to be colored.         | `None`        |
| `colors`             | Sequence or mapping specifying colors of clades.     | `None`        |
| `ax`                 | Axes on which to draw.                               | `None`        |

The least straightforward of these parameters are `clades` and
`colors`. Together, the two parameters allow you to define how certain clades
should be colored. `clades` should map some values identifying clades (clade
names, for example) to the clades they identify.

`colors` can either be a `matplotlib.colors.ListedColormap`, a sequence of RGB
values, or a mapping to RGB values. In the first two cases, the clades will
simply be assigned values from the colormap or the sequence in the order they
appear in the `clades` `dict`. In the last case, each key-value pair specifies
that the clade associated with the same key in the `clades` `dict` should be
colored as specified by the value in the `colors` `dict`.

#### get\_clades

It is sometimes desirable to color clades such that only maximal clades such
that all samples in the clade have some metadata value, *and* all samples that
have that metadata value fall in that clade. To find such clades to be colored,
you can use the `get_clades` function from the `phylo_utils.py` module. 

`get_clades` returns the axes on which the tree was drawn. The parameters
accepted by `get_clades` are described below.

| Formal parameter     | Description                                              |
|----------------------|----------------------------------------------------------|
| `tree`               | The tree in which to find clades.                        |
| `sample_metadata`    | The metadata associated with the samples (nonterminals). |
| `sample_name_column` | Column of `sample_metadata` containing sample names.     |
| `group_by`           | Columns of metadata for which to find clades.            |

```python
import Bio.Phylo
import Bio.Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from phylo_utils import tril_jagged, draw_tree, get_clades
from path_to_sample import path_to_sample

# Convert RNA-clique distance matrix to a 
# Bio.Phylo.TreeConstruction.DistanceMatrix.
phylo_dis_mat = DistanceMatrix(
    [path_to_sample(x) for x in sim.samples],
	tril_jagged(sim.get_dissimilarity_matrix())
)
# Build a tree from the distance matrix using neighbor-joining algorithm.
constructor = DistanceTreeConstructor()
tree = constructor.nj(phylo_dis_mat)
# Get clades by population.
clades = get_clades(tree, sample_metadata, "sample_name", ["population"])
# Draw the tree, coloring maximal clades containingly only and all samples from
# a given population.
ax = draw_tree(
    tree,
	clades=clades,
	colors=matplotlib.colormaps.get("tab10")
)
```

#### draw\_clade\_labels

BioPython can draw clade labels directly on nonterminals of the tree when the
nonterminal is labeled via its `name` attribute, but the author of RNA-clique
has found that labels are often better shown using "calipers"/"braces" on the
side of the plot. To draw clade labels with calipers, you can use the
`draw_clade_labels` function in the `phylo_utils.py` module.

The parameters accepted by `draw_clade_labels` are shown below.

| Formal parameter | Description                                               | Default value |
|------------------|-----------------------------------------------------------|---------------|
| `ax`             | `matplotlib.axes.Axes` on which to draw the labels.       |               |
| `clades`         | Clades for which to draw labels.                          |               |
| `colors`         | Sequence or mapping specifying how to color clade labels. |               |
| `line_padding`   | Horizontal padding for calipers.                          | `0.036`       |
| `cap_width`      | Width of horizontal "caps" on calipers.                   | `0.02`        |
| `text_padding`   | Horizontal padding on text labels next to calipers.       | `0.023`       |
| `make_label`     | Function to make clade label from key in `clades`.        | `id_`         |

`clades` and `colors` work as they do in `draw_tree`. `clades` maps some keys
identifying clades to the clades themselves, and `colors` maps those same keys
to the colors of the clades.

If you don't want the key in the `clades` `dict` to be the label shown on the
plot for that clade, use the `make_label` parameter. By default, the
`make_label` parameter is set to an identity function, so the displayed labels
are exactly the keys in the `clades` `dict`.

```python
from matplotlib import pyplot as plt
from phylo_utils import make_clade_labels

# Make clade labels using the same clades and colors passed to draw_tree.
# Labels should be in all caps.
make_clade_labels(
	ax,
	clades,
	colors,
	make_label=str.upper
1)
```

A real example of using `draw_tree` can also be found at
[`docs/tutorials/reads2tree/make_tree.py`](docs/tutorials/reads2tree/make_tree.py).
