import argparse
import re
import pickle
import pandas as pd
import sys
import itertools
import config as config_module

from pathlib import Path

from tqdm import tqdm

from subset_comparisons import (
    handle_filters,
    matcher,
    make_subset_comparisons,
    relative_to
)
from path_to_sample import sample_re
from build_graph import build_graph
from gene_matches_tables import read_table
from find_homologs import eprint

from collections.abc import Iterable
from typing import Iterator

from gene_matches_tables import get_table_files

def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager()
    arg_config.expose_fields_with_default_aliases(
        "tables_dir",
        "graph",
        required=True
    )
    arg_config.expose_fields_with_default_aliases("output_dir", "title")
    arg_config.expose_config_field("subset_of", aliases=["-I"])
    arg_config.set_required("path_to_sample")
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
    inputs = list(get_table_files(super_config.tables_dir))    
    if args.show_included:
        for sample in super_config.path_to_sample.values():
            if matches(sample):
                print(sample)
        # for df_path in inputs:
        #     df = read_table(df_path, head=1)
        #     if args.show_included and all(
        #         matches(
        #             config.path_to_sample[
        #                 Path(df[x + "sample"].iloc[0]).name
        #             ]
        #         )
        #         for x in ["q", "s"]
        #     ):
        #         print(df_path)
    else:
        #config_module.RNACliqueConfigArgumentManager.make_output_dirs(config)
        if config.tables_dir is not None and \
           config.tables_dir != super_config.tables_dir:
            eprint(
                ("tables_dir is {} but should not be set for this "
                 "program.").format(repr(config.tables_dir))
            )
            eprint("Failing.")
            sys.exit(1)
        config.tables_dir = super_config.tables_dir
        #config.output_dir.mkdir(exist_ok=True)
        config.path_to_sample = {
            p: s for (p, s) in super_config.path_to_sample.items() if matches(p)
        }
        graph = build_graph(
            make_subset_comparisons(
                tqdm(inputs),
                config.tables_dir,
                config.path_to_sample.__contains__
                # matches,
                # config.path_to_sample.__getitem__,
            )
        )
        with open(config.graph / "graph.pkl", "wb") as f:
            pickle.dump(graph, f, pickle.HIGHEST_PROTOCOL)
        config.mark_finish()
        config.yaml_save(args.output_config)


if __name__ == "__main__":
    main()


    
