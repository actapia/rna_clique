import Bio.SeqIO
import config as config_module
from select_top_genes import TopGeneSelector
from transcripts import TranscriptID

from pathlib import Path
from joblib import Parallel, delayed

def select_top_and_save(out_dir, transcripts, x: Path, *args):
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
    arg_config = config_module.RNACliqueConfigArgumentManager()
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
