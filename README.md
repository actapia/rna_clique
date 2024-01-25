# RNA-clique

This is the repository for RNA-clique, a tool for computing pairwise genetic distances from RNA-seq data. The software accepts as input assembled transcriptomes from two or more samples and produces as its output a matrix containing pairwise distances ranging from 0 to 1.

## Installation

This software is written in Python, Perl, and Bash. The software additionally requires NCBI BLAST+ and several Python and Perl libraries. The section below lists the software requirements, and [guides](#installation-guides) are provided for installation of specific systems.

### Requirements

The requirements below represent one tested configuration. This software may work with different versions of the dependencies listed, but such configurations are considered untested.

#### Main software

* Python 3.11
* Perl 5.36.0
* Bash 5.2.15
* ncbi-blast 2.12.0+
* Python libraries
  * tqdm
  * pandas
  * joblib
  * networkx
  * [simple-blast](https://github.com/actapia/simple_blast)
* Perl libraries
  * Bio::SeqIO
  
  
#### Phylogenetics and Visualization

* Python libraries
  * BioPython
  * matplotlib
  * seaborn
  

### Installation guides

* [Ubuntu](
