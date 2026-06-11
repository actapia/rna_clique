import multiprocessing
import sys

from pathlib import Path
from typing import Callable, Iterable, Optional

from multiset_key_dict import MultisetKeyDict

from . import config as config_module
from . import filtering_step
from . import app
from .transcripts import (
    TranscriptID,
    TranscriptIDParseError,
    default_gene_re
)
from .filtered_distance import SampleSimilarity, NoIdealComponentsError
from .similarity_computer import ComparisonSimilarityComputer
from .app import eprint, validate_input_dirs, set_except_hook

def build_parser():
    parser = filtering_step.build_parser()
    parser.parser.description = ("Get a genetic distance matrix from input "
                                 "transcriptomes.")
    parser.expose_fields_with_default_aliases("matrix")
    return parser

def rna_clique(
        dirs: Iterable[Path],
        out_dir_1: Path,
        out_dir_2: Path,
        cache_dir: Path,
        output_graph: Path,
        output_matrix: Optional[Path],
        top_genes: int,
        transcripts: str = "transcripts.fasta",
        top_matches: int = 1,
        id_parser: Callable[
            [str],
            TranscriptID
        ] = TranscriptID.parser_from_re(default_gene_re),
        evalue: float = 1e-99,
        keep_all: bool = True,
        store_dfs: bool = False,
        jobs: int = multiprocessing.cpu_count() - 1,
) -> tuple[SampleSimilarity, dict[Path, str]]:
    """Perform a full RNA-clique analysis using the provided transcriptomes.

    This function provides an API for performing a full RNA-clique
    analysis---from the first step of selecting top genes by k-mer coverate to
    computing distances from filtered gene matches tables.

    Specifically, the steps this function takes are:

    1. Selecting transcripts of top n genes per sample by k-mer coverage.
    2. Pairwise BLASTing every sample's top n genes against every other's.
    3. Making gene matches tables for each sample pair from BLASTn results.
    4. Building the gene matches graph from the gene matches tables.
    5. Filtering the gene matches tables using the gene matches graph.
    6. Computing similarities (or distances) using filtered tables.

    These steps can also be performed separately using other function available
    in the RNA-clique package.

    This function requires an Iterable of Paths to directories containing
    assembled transcriptomes for the samples to be analyzed. The function also
    requires two output directories---one for the top n genes, and one of the
    gene matches tables---an intermediate directory for storing BLAST DBs, a
    Path at which to store the output gene matches graph, and a path at which to
    store the output distance matrix. Finally, the caller must also provide the
    number of top genes (parameter little n) to select.

    This function also accepts several optional parameters. RNA-clique expects
    all transcriptomes to have the same file name; it distinguishes them by
    their directory names. By default, RNA-clique expects to find transcriptomes
    in files named "tarnscripts.fasta," but the name of this file can be
    changed. The name must, however, be the same for all transcriptomes used.

    RNA-clique produces a given gene matches table for a pair of samples A and B
    using the BLASTn results for comparing A against B and B against A both. For
    a given pair of genes a and b, the match only appears in the resulting gene
    matches table for A and B if an alignment for a against b appears among the
    hits with top N (big N) bitscores for alignments involving a in both
    directions. Big N should usually be set to 1 but is configurable. The effect
    of setting N > 1 has not been carefully tested; it might find more orthologs
    while also introducing some false positive matches that could inflate the
    distance.

    Additionally, when there are more than N pairs of genes such that each pair
    is among the pairs of genes with the N highest bitscores in both directions
    (i.e., when there are ties), the default behavior is to keep all of the gene
    pairs in the resulting gene matches table. Hence, by default, the gene
    matches tables can map a sample A gene to more than N sample B genes. To
    change this behavior, the keep_all parameter can be set to False. In that
    case, ties will be broken arbitrarily, but only N sample B genes will be
    matched to each sample A gene in the gene matches table for A and B.

    RNA-clique expects that each transcript's FASTA ID contains the transcript's
    k-mer coverage, gene ID, and isoform ID. RNA-clique requires a function to
    parse these transcript IDs; this can be provided using the id_parser
    parameter. The provided id_parser should be a function accepting a string
    and returning a TranscriptID object. By default, the parameter uses a
    regex-based parser designed to parse transcript IDs produced by the
    rnaSPAdes assembler. To create a parser for a custom regex, you can use the
    TranscriptID.parser_from_re classmethod.

    By default, the BLASTn searches will use an e-value cutoff of 1e-99; this
    all but guarantees that alignments represent some kind of homology. For less
    closely related samples, it might help to increase the threshold; this can
    be controlled using the evalue parameter.

    RNA-clique ordinarily does not store the gene matches tables in RAM, instead
    opting to use iterators to reduce the memory footprint. If tables are needed
    immediately after running this method, the store_dfs option can be set to
    True to keep the gene matches table dataframes in the SampleSimilarity
    object.

    See the original RNA-clique publication, "RNA-clique: a method for computing
    genetic distances from RNA-seq data" for more details on RNA-clique's
    workings.

    When a non-None value is provided for output_matrix, this function eagerly
    computes the distance matrix and saves it to the provided Path. If there are
    no ideal components in the gene matches graph, computing distances is not
    possible, and attempting to get the distances will raise a
    NoIdealComponentsError. To avoid this possibility, provide None for the
    output_matrix parameter instead.

    This function returns two values. The first is a SampleSimilarity object for
    the analysis. The SampleSimilarity object provides access to the gene
    matches graph, gene matches tables (if stored), and distance matrix. The
    second value returned is a dict mapping the FASTA file Paths for the top n
    genes of each sample to its corresponding sample name.

    Parameters:
        dirs:              Paths to input directories containing transcriptomes.
        out_dir_1:         Path to output directory for top n genes by coverage.
        out_dir_2:         Path to output directory for gene matches tables.
        cache_dir:         Path to BLAST DB intermediate directory.
        output_graph:      Path to output graph.
        output_matrix:     Path to output distance matrix.
        top_genes (int):   Number of top genes by k-mer coverage to select.
        transcripts (str): Name of all transcriptome FASTA files.
        top_matches (int): Order of bitscore needed bidirectionally for a match.
        id_parser:         Function to parse FASTA IDs into TranscriptIDs.
        evalue (float):    e-value threshold to use for BLASTn searches.
        keep_all (bool):   Keep all gene pairs in the case of ties by bitscore.
        store_dfs (bool):  Store gene matches tables in SampleSimilarity object.
        jobs (int):        Number of parallel jobs to use.

    Returns:
        SampleSimilarity with distances and graph and Path-to-sample mapping.
    """
    tables, table_paths, graph, num_tables, pts = filtering_step.filtering_step(
        dirs,
        out_dir_1,
        out_dir_2,
        cache_dir,
        output_graph,
        top_genes,
        transcripts,
        top_matches,
        id_parser,
        evalue,
        keep_all,
        jobs,
    )
    tables = ComparisonSimilarityComputer.mapping_from_dfs(tables)
    if store_dfs:
        tables = MultisetKeyDict(tables)
    sim = SampleSimilarity(
        graph,
        tables,
    )
    if output_matrix is not None:
        mat = sim.get_dissimilarity_df()
        mat.to_hdf(output_matrix, key="matrix")
    return sim, pts
    
def main():
    with set_except_hook():
        _, args, config = build_parser().get_arguments_and_config()
    with set_except_hook(config.verbose):
        config_module.RNACliqueConfigArgumentManager.make_output_dirs(config)
        id_parser = TranscriptID.parser_from_re(config.transcript_id_regex)
        validate_input_dirs(config)
        try:
            sim, pts = rna_clique(
                config.input_dirs,
                config.top_genes_dir,
                config.tables_dir,
                config.cache_dir,
                config.graph,
                None,
                config.top_genes,
                config.transcripts_name,
                config.top_matches,
                id_parser,
                config.evalue,
                config.keep_all,
                False,
                jobs=config.jobs
            )
            config.path_to_sample = pts    
            mat = sim.get_dissimilarity_df()
            mat.to_hdf(config.matrix, key="matrix")
            config.mark_finish()
            if args.output_config:
                config.yaml_save(args.output_config)        
        except NoIdealComponentsError:
            eprint("No ideal components found. Cannot report distances!")
            sys.exit(1)
        except TranscriptIDParseError:
            app.print_transcript_id_parse_error_message(
                config.transcript_id_regex
            )
            raise

    
if __name__ == "__main__":
    main()
