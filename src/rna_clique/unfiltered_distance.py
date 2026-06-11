from . import config as config_module
from .similarity_computer import (
    ComparisonSimilarityComputer,
    similarities_from_dfs
)
from .gene_matches_tables import get_table_files
from .app import set_except_hook, eprint

class UnfilteredSimilarity(ComparisonSimilarityComputer):
    """Computes similarities using gene matches tables, without filtering.
    
    Attributes:
        comparison_dfs: An iterable mapping sample pairs to comparisons.
    """
    def _similarity_helper(self):
        return similarities_from_dfs(self.comparison_dfs)
    
def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager(
        description=(
            "Compute pairwise distasnces from gene matches tables alone."
        ),
    )
    arg_config.expose_fields_with_default_aliases(
        "tables_dir",
        "matrix",    
        required=True
    )
    arg_config.expose_fields_with_default_aliases(
        "output_dir",
    )
    arg_config.add_output_config_argument()
    return arg_config

def main():
    with set_except_hook():
        _, args, config = build_parser().get_arguments_and_config()
    with set_except_hook(config.verbose):
        tables = list(get_table_files(config.tables_dir))
        if not tables:
            eprint(
                "Warning: No gene matches tables found in {}".format(
                    config.tables_dir
                )
            )        
        sim = UnfilteredSimilarity.from_filenames(
            get_table_files(config.tables_dir)
        )
        mat = sim.get_dissimilarity_df()
        mat.to_hdf(config.matrix, key="matrix", mode="w")
        config.mark_finish()
        if args.output_config is not None:
            config.yaml_save(args.output_config)    

if __name__ == "__main__":
    main()

