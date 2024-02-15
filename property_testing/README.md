# Property testing script

The script in this directory tests that distances do not violate symmetry or the
triangle inequality for the given FASTA files. The script is not part of the
main RNA-clique program but was used for results reported in the paper
*RNA-clique: A method for computing genetic distances from RNA-seq data*.

## Usage

This script expects to be in the root of the RNA-clique repository. Before using
it, you can make a link to the script from the root.

```bash
ln -s property_testing/test_properties.py .
```
