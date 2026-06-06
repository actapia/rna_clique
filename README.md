# RNA-clique

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14890599.svg)](https://doi.org/10.5281/zenodo.14890599)

This is the repository for RNA-clique, a tool for computing pairwise genetic
distances from RNA-seq data. The software accepts as input assembled
transcriptomes from two or more samples and produces as its output a matrix
containing pairwise distances ranging from 0 to 1.

## Installation

This software is written in Python. The software additionally requires NCBI
BLAST+ and several Python libraries. [Guides](#installation-guides) are provided
for installation on specific systems. Alternatively, for installing on other
systems, you can see the
<!--{{doc_link("requirements.md", "requirements", False) | comment_surround}}{{empty("-->
[requirements](https://actapia.github.io/rna_clique/dev/requirements)<!--")}}-->

### Installation guides

<!--{{doc_link("installation_guides/ubuntu.md", "Ubuntu", False) | comment_surround}}{{empty("-->
* [Ubuntu](https://actapia.github.io/rna_clique/dev/installation_guides/ubuntu)
<!--")}}-->
<!--{{doc_link("installation_guides/macos.md", "macOS", False) | comment_surround}}{{empty("-->
* [macOS](https://actapia.github.io/rna_clique/dev/installation_guides/macos)
<!--")}}-->

### Basic usage

To run RNA-clique on your assembled transcriptomes, first make sure that your
data are in a <!--{{doc_link("formats.md#transcriptomes", "format understood by
RNA-clique", False) | comment_surround}}{{empty("-->
[format understood by RNA-clique](https://actapia.github.io/rna_clique/dev/formats.md#transcriptomes).
<!--")}}-->

Then, run `rna_clique.py` from the root of this repository with the directories
containing your transcriptomes, an output directory, and a setting for the
number of top genes to select.

```bash
python rna_clique.py -O my_rna_clique_out -n 50000 \
                     path/to/transcriptome1_dir \
                     path/to/transcriptome2_dir \
					 path/to/transcriptome3_dir ...
```

RNA-clique produces an output matrix at `my_rna_clique_out/matrix.h5`. To see it
in a human-readable format, use `export_matrix.py`.

```bash
python export_matrix -m my_rna_clique_out/matrix.h5 
```

More details about the usage of RNA-clique can be found in the <!--{{doc_link("usage.md", "Command-line usage guide", False) | comment_surround}}{{empty("-->[Command-line usage guide](https://actapia.github.io/rna_clique/dev/usage)<!--")}}-->


### Downstream analyses

The `export_matrix.py` script prints the calculated matrix to the standard
output, so you can use redirection or pipes to save the results to a file. You
could then use the matrix in any downstream application capable of loading
arbitrary matrices from files.

For example, if you output the matrix to a file named `distances`, you could
load the matrix in R using the following code:

```R
dis <- as.matrix(read.table("distances", sep=" "))
```

## Using RNA-clique in Python code

You can use RNA-clique directly from your Python code. For example,

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

For information on finer-grained control via RNA-clique's Python API, see the
<!--{{doc_link("api/README.md", "API Guide", False) | comment_surround}}{{empty("-->
[API guide](https://actapia.github.io/rna_clique/dev/api).
<!--")}}-->

## License

All code except tests/test_install/macos.sh is licensed under the MIT license,
which may be found at LICENSE.MIT at the root of this repository.

[tests/test_install/macos.sh](https://github.com/actapia/rna_clique/blob/main/tests/test_install/macos.sh)
is licensed under the Creative Commons Attribution-ShareAlike 4.0 License, which
may be found at LICENSE.CC-BY-SA-4.0 at the root of this repository.

A machine-readable copyright file in Debian format may also be found at
[copyright](https://github.com/actapia/rna_clique/blob/main/copyright).

## Citation

If you use RNA-clique for your work, please cite ["RNA-clique: a method for
computing genetic distances from RNA-seq
data"](https://doi.org/10.1186/s12859-024-05811-9).

<!-- {% raw %}{{ -->
```tex
@article{tapia2024rna,
  title={{RNA-clique: a method for computing genetic distances from RNA-seq data}},
  author={Tapia, Andrew C and Jaromczyk, Jerzy W and Moore, Neil and Schardl, Christopher L},
  journal={BMC Bioinformatics},
  volume={25},
  year={2024},
  publisher={BioMed Central},
  keywords={pub}
}
```
<!-- }}{% endraw %} -->

## Additional documentation

* <!--{{doc_link("usage.md", "Command-line usage guide", False) | comment_surround}}{{empty("-->[Command-line usage guide](https://actapia.github.io/rna_clique/dev/usage)<!--")}}-->
* <!--{{doc_link("tutorials/reads2tree/README.md", "Tutorial: From RNA-seq reads to a phylogenetic tree with RNA-clique", False) | comment_surround}}{{empty("-->[Tutorial: From RNA-seq reads to a phylogenetic tree with RNA-clique](https://actapia.github.io/rna_clique/dev/tutorials/reads2tree)<!--")}}-->
