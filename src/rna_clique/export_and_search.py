import argparse
import json
import sys
import re
import app

import search_ideal_components
import export_orthologs
import config as config_module

from name_conflict_resolver import NameConflictResolver

from collections import Counter
from collections.abc import Iterable
from typing import Optional, Callable
from pathlib import Path

from tqdm import tqdm

from filtered_distance import SampleSimilarity, NoIdealComponentsError
from gene_matches_tables import get_table_files
from transcripts import default_gene_re, TranscriptID, TranscriptIDParseError
from app import set_except_hook, eprint
from path_to_sample import (
    path_to_sample,
    dict_path_to_sample,
    sample_re,
    PathToSampleError
)

# @marshalling_dataclass(optional=True)
# class ExportAndSearchConfig:
#     """Dataclass with config data for the export and search program.."""
#     config_version: str = marshalling_field(default="0.0.1", metadata={
#         "description": "Version of the configuration schema used."})
#     resolve_name_conflicts: Optional[bool] = marshalling_field(metadata={
#         "description": "Automatically resolve conflicting filenames."
#     })
#     export_output_dir: Optional[Path] = marshalling_field(metadata={
#         "description": "Output directory for exported sequences."
#     })    
    

def build_parser():
    parser = config_module.ArgumentManager(
        description=(
            "Export orthologs from ideal components and search their sequences."
        ),
    )
    #parser.add_argument("-A", "--analyses", type=Path, nargs="+", required=True)
    parser.add_argument(
        "-C",
        "--configs",
        type=Path,
        nargs="+",
        required=True,
        help="RNA-clique configs for which to export and search orthologs.",
    )
    parser.add_argument(
        "-r",
        "--resolve-name-conflicts",
        action="store_true",
        help="Resolve conflicting output filenames automatically.",
    )
    parser.add_argument(
        "-X",
        "--export-output-dir",
        type=Path,
        required=True,
        help="Output directory for exported sequences.",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        help="Number of parallel jobs to use.",
    )
    parser.add_argument(
        "-x",
        "--export-only",
        action="store_true",
        help="Only export the orthologs; don't search."
    )
    parser.add_argument(
        "-Q",
        "--queries",
        nargs="+",
        type=Path,
        help="FASTA files containing sequences to search in orthologs.",
    )
    parser.add_argument(
        "--transcript-id-regex",
        "-p",
        "--pattern",
        type=re.compile,
        help="Python regex for parsing sequence IDs",
        default=default_gene_re
    )
    parser.add_argument(
        "--extended-search-evalue",
        "-E",
        const=search_ideal_components.default_extended_search_evalue,
        nargs="?",
        type=float,
        help="Search other isoforms of a gene that produces a hit.",
    )
    parser.add_argument(
        "--search-evalue",
        "-e",
        type=float,
        help="e-value cutoff to use for initial searches.",
        default=search_ideal_components.default_search_evalue,
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Print more output than ususal."
    )
    return parser

def handle_arguments():    
    return build_parser().parse_args()

def get_analysis_name(config: config_module.RNACliqueConfig) -> str:
    """Get a suitable name for the analysis represented by the configuration.

    If the configuration includes a title, this is preferred. Otherwise, the
    name of the output directory is used. If neither is present in the
    configuration, this function raises a ValueError.

    Parasmeters:
        config: RNACliqueConfig representing analysis for which to get a name.

    Returns:
        A string that can be used as a name for the analysis.
    """
    name = config.title
    if not name:
        try:
            out_dir = config.output_dir
        except AttributeError:
            raise ValueError("Could not determine analysis name.")
        name = out_dir.name
    return name

class NameConflictError(ValueError):
    pass

class RegexTranscriptIDParseError(TranscriptIDParseError):
    def __init__(self, message, regex):
        super().__init__(message)
        self.regex = regex

def export_and_search(
        configs: list[config_module.RNACliqueConfig],
        export_output_dir: Path,
        queries: Iterable[Path],
        parse_transcript_id: Optional[Callable[[str], TranscriptID]] = None,
        jobs: Optional[int] = None,
        resolve_name_conflicts: bool = False,
        export_only: bool = False,
        extended_evalue: Optional[bool | float] = None,
        evalue: float = search_ideal_components.default_search_evalue,
):
    """Export transcripts in ideal components and search them for sequences.

    This function allows you to export transcripts belonging to genes in ideal
    component and search for sequences in multiple query FASTA files within
    them for multiple analyses at once. Each analysis should be represented by
    an RNACliqueConfig object.

    Parameters:
        configs (list):                Configs for analyses to export/search.
        export_output_dir:             Export/search result output directory.
        queries:                       Paths to query sequence FASTA files.
        parse_transcript_id:           Function to parse transcript FASTA IDs.
        jobs (int):                    Number of parallel jobs to use.
        resolve_name_conflicts (bool): Auto-resolve conflicting analysis names.
        export_only (bool):            Only perform the export step.
        extended_evalue (float):       e-value cutoff for extended searches.
        evalue (float):                e-value cutoff for initial searches.
    """
    out_names = [export_output_dir / get_analysis_name(c) for c in configs]
    counts = Counter(out_names)
    try:
        multiple = next(x for (x, c) in counts.items() if c > 1).name
        #eprint(f"Multiple analyses named {multiple!r}.")
        if not resolve_name_conflicts:
            raise NameConflictError(f"Multiple analyses named {multiple!r}.")
        # else:
        #     eprint("Resolving conflicts automatically.")
    except StopIteration:
        pass
    resolver = NameConflictResolver.from_keys(
        range(len(out_names)),
        out_names.__getitem__
    )
    for config, (_, (_, out_dir)) in zip(configs, resolver.resolve()):
        # if config.transcript_id_regex is None:
        #     config.transcript_id_regex = transcript_id_regex
        if jobs is not None:
            config.jobs = jobs
        out_dir.mkdir(exist_ok=True)
        export_dir = out_dir / "export"
        export_dir.mkdir(exist_ok=True)
        comparison_paths = list(get_table_files(config.tables_dir))
        sim = SampleSimilarity.from_filenames(
            config.graph,
            tqdm(comparison_paths),
            store_dfs=True
        )
        try:
            sim.similarities
        except NoIdealComponentsError:
            eprint(
                "Warning: One or more analyses has no ideal components to "
                "export! Skipping ... "
            )
            continue
        pti = parse_transcript_id
        if pti is None:
            pti = TranscriptID.parser_from_re(
                config.transcript_id_regex,
            )
        if config.path_to_sample is not None:
            pts = dict_path_to_sample(config.path_to_sample)
        else:
            pts = path_to_sample
        exporter = export_orthologs.OrthologExporter(
            sim,
            pti,
            False,
            debug=True,
            consistent_strands=True,
            allow_inconsistent=True,
            jobs=config.jobs,
            path_to_sample=pts,
        )
        component_paths = exporter.by_component(export_dir, order="after")
        if not export_only:
            all_ideal_path = export_dir / "all_ideal.fasta"
            # with open(all_ideal_path, "w") as all_ideal:
            #     for path in component_paths.values():
            #         with open(path, "r") as component_fasta:
            #             all_ideal.write(component_fasta.read())
            exporter.make_all_ideal(component_paths, export_dir)
            db_cache = export_dir / "db_cache"
            db_cache.mkdir(exist_ok=True)
            for query in queries:
                search_dir = out_dir / ("search_" + query.stem)
                search_dir.mkdir(exist_ok=True)
                stats = search_ideal_components.search(
                    sim=sim,                    
                    exported=all_ideal_path,
                    db_cache_loc=db_cache,
                    out_dir=search_dir,
                    query=query,
                    merge_sams=True,
                    parse_transcript_id=pti,
                    jobs=config.jobs,
                    strand_graph=exporter.strand_graph,
                    node_to_ccc=exporter.node_to_component_component,
                    extended_evalue=extended_evalue,
                    evalue=evalue,
                    path_to_sample=pts
                )
                # if stats is None:
                #     from IPython import embed; embed()
                with open(search_dir / "stats", "w") as stats_file:
                    json.dump(stats._asdict(), stats_file)

def main():
    with set_except_hook():
        _, args = build_parser().get_arguments()
    with set_except_hook(args.verbose):
        for query in args.queries:
            config_module.RNACliqueConfig.validate_file(query, "Query file ")
            search_ideal_components.check_is_not_empty(query, "query")
        configs = [
            config_module.RNACliqueConfig.yaml_load(c) for c in args.configs
        ]
        no_path_to_sample = [
            p for (c, p) in zip(configs, args.configs)
            if c.path_to_sample is None
        ]
        if no_path_to_sample:
            eprint("Warning: the following configuration files lack "
                   "path_to_sample attributes:\n")
            for p in no_path_to_sample:
                eprint(p)
            eprint("\nFor those analyses, this script will parse sample names "
                   "from paths using the default regular expression "
                   f"{sample_re.pattern}.")                
        args.export_output_dir.mkdir(exist_ok=True)
        parse_transcript_id = TranscriptID.parser_from_re(
            args.transcript_id_regex
        )
        # if .path_to_sample is not None:
        #     pts = dict_path_to_sample(config.path_to_sample)
        #     mode = f"path_to_sample in configuration file {args.input_config}"
        # else:
        #     eprint("Warning: path_to_sample attribute not found in input "
        #            "config. Will parse sample names from paths using default "
        #            f"regular expression {sample_re.pattern}.")
        #     pts = path_to_sample
        #     mode = f"sample_regex {sample_re.pattern}"        
        try:
            export_and_search(
                configs,
                args.export_output_dir,
                args.queries,
                parse_transcript_id,
                args.jobs,
                args.resolve_name_conflicts,
                args.export_only,
                args.extended_search_evalue,
                args.search_evalue,
            )
        except NameConflictError as e:
            eprint(e)
            eprint(("Cannot continue. Provide the -r option to try "
                    "automatic resolution."))
            sys.exit(1)
        except PathToSampleError as e:
            eprint(
                f"RNA-clique could not get the sample name for path {e.path}. "
                "Please ensure that path_to_sample is provided in the configs "
                "for any analyses where sample names cannot be parsed with the "
                f"default regex {sample_re.pattern}."
            )
            raise e
        except TranscriptIDParseError as e:
            app.print_transcript_id_parse_error_message(
                args.transcript_id_regex
            )
            raise e
        except export_orthologs.ExportTooManyFilesError as e:
            app.print_too_many_files_error_message(e)
            raise e
        
if __name__ == "__main__":
    main()
