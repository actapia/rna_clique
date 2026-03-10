import re
import pickle
import sys
import copy

import config as config_module

from pathlib import Path
from typing import Callable

from tqdm import tqdm

from subset_comparisons import (
    handle_filters,
    matcher,
    make_subset_comparisons,
)
from build_graph import build_graph
from find_homologs import eprint
from gene_matches_tables import get_table_files

def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager(
        description=(
            "Create an analysis by reusing previous results for a subset of "
            "samples."
        ),
    )
    arg_config.expose_fields_with_default_aliases(
        "tables_dir",
        "graph",
        required=True
    )
    arg_config.expose_fields_with_default_aliases("output_dir", "title")
    arg_config.expose_config_field("subset_of", aliases=["-I"], required=True)
    arg_config.set_defaults("top_genes_dir", None)
    #arg_config.set_required("path_to_sample")
    arg_config.add_argument(
        "--exclude",
        "-x",
        nargs="+",
        default=[],
        help="samples to exclude (default is none)"
    )
    arg_config.add_argument(
        "--include",
        "-y",
        nargs="+",
        default=[],
        help="samples to include (default is all)"
    )
    arg_config.add_argument(
        "--include-regex",
        "-Y",
        type=re.compile,
        help="regular expression specifying which sample names to include"
    )
    arg_config.add_output_config_argument()
    # arg_config.add_argument(
    #     "--sample-name-regex",
    #     "-r",
    #     help="regular expression to apply to sample names",
    #     type=re.compile,
    #     default=sample_re
    # )
    arg_config.add_argument(
        "--include-file",
        type=Path,
        help="file containing samples to include"
    )
    arg_config.add_argument(
        "--exclude-file",
        type=Path,
        help="file containing samples to exclude"
    )
    arg_config.add_argument(
        "--show-included",
        action="store_true",
        help="show which samples would be included and exit"
    )
    arg_config.set_defaults("top_genes_dir", None)
    # arg_config.add_argument(
    #     "--show-parsed-paths",
    #     action="store_true",
    #     help="show parsed paths"
    # )
    return arg_config

class SubsetAnalysisCreator:
    """Class for creating analysis whose samples are a subset of another's.

    For any RNA-clique anaylsis A with a set of S samples and gene matches
    tables G and analysis B with a set of T samples and gene matches tables H,
    if T is a subset of S, then H must be a subset of G. This means that we can
    avoid recomputing gene matches tables if we wish to obtain distances for a
    subset of an existing analysis, and this can save much time. This class
    facilitates creating such subset analyses.

    Much of this class performs I/O operations, but the class also provides a
    config property that provides an RNACliqueConfig object for the new
    analysis.

    Attributes:
        matches:      Predicate indicating if sample is included in subset.
        super_config: Configuration for the parent analysis to be subsetted.
        config:       Configuration for the child (subset) analysis.
    """
    def __init__(
            self,
            matches: Callable[[Path], bool],
            super_config: config_module.RNACliqueConfig,
            config: config_module.RNACliqueConfig,
    ):
        """Construct a SubsetAnalysisCreator with predicate and configs.

        The first argument, matches, should be a predicate on Paths indicating
        whether its argument is to be included in the subset.        

        Note that a copy of the provided config is made and set for the
        created SubsetAnalysisCreator.

        It is expected that the child configuration contains non-None tables_dir
        and graph attributes. If these are not found, the constructor raises a
        ValueError.

        Parameters:
            matches:      Predicate indicating if sample should be included.
            super_config: Configuration of parent analysis.
            config:       Configuration of the child analysis.
        """
        if config.tables_dir is None:
            raise ValueError("tables_dir attribute of child config is None.")
        if config.graph is None:
            raise ValueError("graph attribute of child config is None.")        
        self.matches = matches
        self.super_config = super_config
        # self.tables_dir = tables_dir
        # self.graph = graph
        self.config = copy.deepcopy(config)

    @classmethod
    def from_paths(
            cls,
            matches: Callable[[Path], bool],
            super_config: config_module.RNACliqueConfig,
            tables_dir: Path,
            graph: Path,            
    ):
        """Construct a SubsetAnalysisCreator from output paths.

        This classmethod allows one to construct a SubsetAnalysisCreator without
        explicitly constructing a new RNACliqueConfig first.

        Parameters:
            matches:      Predicate indicating if sample should be included.
            super_config: Configuration of parent analysis.
            tables_dir:   Output path for gene matches table links.
            graph:        Output path for new gene matches graph.
        """
        return cls(
            matches,
            super_config,
            config_module.RNACliqueConfig(tables_dir=tables_dir, graph=graph)
        )

    def make(self):
        """Make the child (subset) analysis, linking tables and making a graph.

        This method creates symlinks to the parent analysis's gene matches
        tables within the directory specified by the tables_dir attribute of the
        child config. It then builds a graph using the subset of tables.

        This method also updates the SubsetAnalysisCreator's config attribute
        using values from the parent config.
        """
        inputs = list(get_table_files(self.super_config.tables_dir))
        self.config.top_genes_dir = self.super_config.top_genes_dir
        self.config.tables_dir.mkdir(exist_ok=True)
        self.config.path_to_sample = {
            p: s for (p, s) in self.super_config.path_to_sample.items()
            if self.matches(s)
        }
        self.config.input_dirs = [
            p for p in self.super_config.input_dirs if self.matches(p.name)
        ]
        self.config.top_genes = self.super_config.top_genes
        self.config.transcripts_name = self.super_config.transcripts_name
        self.config.top_matches = self.super_config.top_matches
        self.config.evalue = self.super_config.evalue
        self.config.keep_all = self.super_config.keep_all
        self.config.jobs = self.super_config.jobs
        self.config.transcript_id_regex = self.super_config.transcript_id_regex
        graph = build_graph(
            make_subset_comparisons(
                tqdm(inputs),
                self.config.tables_dir,
                self.config.path_to_sample.__contains__
            )
        )
        with open(self.config.graph, "wb") as f:
            pickle.dump(graph, f, pickle.HIGHEST_PROTOCOL)        

def main():
    _, args, config = build_parser().get_arguments_and_config()
    include = handle_filters(args.include, args.include_file)
    if not include:
        include = None
    exclude = handle_filters(args.exclude, args.exclude_file)
    if not exclude:
        exclude = None
    matches = matcher(
        include,
        exclude,
        args.include_regex
    )
    super_config = config_module.RNACliqueConfig.yaml_load(args.subset_of)
    if super_config.path_to_sample is None:
        eprint("Parent config must contain path_to_sample attribute.")
        sys.exit(1)
    if args.show_included:
        for sample in super_config.path_to_sample.values():
            if matches(sample):
                print(sample)
    else:
        if config.top_genes_dir is not None and \
           config.top_genes_dir != super_config.top_genes_dir:
            eprint(
                ("top_genes_dir is {} but should not be set for this "
                 "program.").format(repr(config.top_genes__dir))
            )
            eprint("Failing.")
            sys.exit(1)
        if config.output_dir is not None:
            config.output_dir.mkdir(exist_ok=True)
        creator = SubsetAnalysisCreator(matcher, super_config, config)
        creator.make()
        creator.config.mark_finish()
        creator.config.yaml_save(args.output_config)

if __name__ == "__main__":
    main()


    
