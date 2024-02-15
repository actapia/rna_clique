# Parameter testing scripts

The scripts in this directory were used for the parameter tests described in
*RNA-clique: A method for computing genetic distances from RNA-seq data* (plus
some tests that are not described there) but are not part of the main RNA-clique
program. Hence, detailed usage information about these script is not (yet)
provided here.

## Usage

Most scripts in this directory expect to be in the root of the RNA-clique
directory. Before using them, run a command like the following from the root of
the repository.

```bash
for f in parameter_testing/*.{py,sh}; do ln -s "$f" .; done
```


## Files

| Name                                    | Description                                                                                                                                                          |
|-----------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `add_output_params.py`                  | Processes command-line arguments for `stats_for_params.sh` to add output directories automatically based on existing command-line arguments.                         |
| `do_sample_count_tests.sh`              | Performs tests assessing the effect of sample count ($s$) on component counts for prefixes of a randomly selected permutation of the input FASTA files.              |
| `do_top_n_tests.sh`                     | Performs tests assessing the effects of top genes selected ($n$) on components counts for the specified input FASTA files.                                           |
| `fair_sample_count_tests.py`            | Performs tests assessing the effect of sample count ($s$) and top genes selected ($n$) on component counts for many random subsets of the set of input FASTA files.  |
| `gen_params.py`                         | Generates command-line arguments for `stats_for_params.sh`. Arguments to `gen_params.py` specify which options to generate. Lines of standard input specify  values. |
| `genotype_endophyte_interleave_test.py` | Performs tests assessing the effect of $s$ on component counts for prefixes of a permutation of the input FASTA file interleaved by genotype + endophyte status.     |
| `genotype_interleave_test.py`           | Performs tests assessing the effect of sample count ($s$) on component counts for prefixes of a permutation of the input FASTA files interleaved by genotype.        |
| `genotype_order_test.py`                | Performs tests assessing the effect of sample count ($s$) on component counts for prefixes of a permutation of the input FASTA files ordered by genotype.            |
| `incseq.sh`                             | Outputs an arithmetic sequence, including endpoints.                                                                                                                 |
| `one_nonzero.py`                        | Outputs a sequence based on the OEIS A037124 "one nonzero" sequence, i.e., 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40 ...                                             |
| `ordered_sample_count_test.py`          | Python module containing a class for permutation prefix sample count ($s$) tests.                                                                                    |
| `sample_count_viz_utils.py`             | Utilities used for visualizing sample count test results.                                                                                                            |
| `stats_for_params.sh`                   | Performs phase 1 of RNA-clique with the provided parameters and computes component count statistics.                                                                 |
| `try_params.sh`                         | Runs `stats_for_params.sh` for various parameters settings. Options to be set are provided as arguments. Values are provided as lines in standard input.             |
