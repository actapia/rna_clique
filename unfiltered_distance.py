import config as config_module

from similarity_computer import (
    ComparisonSimilarityComputer,
    similarities_from_dfs
)
from gene_matches_tables import get_table_files

class UnfilteredSimilarity(ComparisonSimilarityComputer):
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
    _, args, config = build_parser().get_arguments_and_config()
    sim = UnfilteredSimilarity.from_filenames(
        get_table_files(config.tables_dir)
    )
    mat = sim.get_dissimilarity_df()
    mat.to_hdf(config.matrix, key="matrix", mode="w")
    config.mark_finish()
    config.yaml_save(args.output_config)    

if __name__ == "__main__":
    main()

