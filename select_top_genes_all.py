import Bio.SeqIO

import config as config_module

from pathlib import Path
from joblib import Parallel, delayed

from select_top_genes import TopGeneSelector
from transcripts import TranscriptID

def select_top_and_save(
        out_dir: Path,
        transcripts: str,
        x: Path,
        *args
) -> tuple[Path, str]:
    """Select the top n genes by k-mer coverage from transcripts and save them.

    This function wraps TopGeneSelector, saving the selected top transcripts by
    k-mer coverage in a specified directory. The name given to the output file
    is "{sample}_top.fasta", where the sample is assumed to be the name of the
    directory in which the transcripts FASTA file is located.

    Additional arguments provided will be passed to the
    TopGeneSelector.from_path classmethod used to construct a TopGeneSelector
    object.

    Parameters:
        out_dir:           Location in which to save top n genes.
        transcripts (str): Name of the FASTA file containing transcripts.
        x:                 Directory containing the transcripts FASTA file.

    Returns:
        Path to the output file and the inferred sample name.
    """
    out = out_dir / (x.stem + "_top.fasta")
    Bio.SeqIO.write(
        TopGeneSelector.from_path(
            x / transcripts,
            *args
        ).get_top_gene_seqs(),
        out,
        "fasta"
    )
    return (out, x.stem)

def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager(
        description="Select top genes by k-mer coverage for given samples.",
    )
    arg_config.expose_fields_with_default_aliases(
        "top_genes",
        "transcripts_name",
        "top_genes_dir",
        "transcript_id_regex",
        "jobs",
        required=True
    )
    arg_config.expose_fields_with_default_aliases("output_dir", "title")
    arg_config.add_output_config_argument()
    arg_config.expose_config_field(
        "input_dirs",
        nargs="*",
        #const=None,
        positional=True,
    )
    return arg_config

def main():
    _, args, config = build_parser().get_arguments_and_config()
    config_module.RNACliqueConfigArgumentManager.make_output_dirs(config)
    id_parser = TranscriptID.parser_from_re(config.transcript_id_regex)
    pts = dict(
        Parallel(n_jobs=config.jobs)(
            map(
                delayed(
                    lambda x: select_top_and_save(
                        config.top_genes_dir,
                        config.transcripts_name,
                        x,
                        config.top_genes,
                        id_parser
                    )
                ),
                args.input_dirs
            )
        )
    )
    config.path_to_sample = pts
    config.mark_finish()
    config.yaml_save(args.output_config)

if __name__ == "__main__":
    main()
