import argparse
import json
import sys
import re

import search_ideal_components
import export_orthologs
import config as config_module

from name_conflict_resolver import NameConflictResolver

from collections import Counter
from collections.abc import Iterable
from typing import Optional, Callable
from pathlib import Path

from tqdm import tqdm

from filtered_distance import SampleSimilarity
from find_homologs import eprint
from gene_matches_tables import get_table_files
from transcripts import default_gene_re, TranscriptID

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

def export_and_search(
        configs: list[config_module.RNACliqueConfig],
        export_output_dir: Path,
        queries: Iterable[Path],
        parse_transcript_id: Optional[Callable[[str], TranscriptID]] = None,
        jobs: Optional[int] = None,
        resolve_name_conflicts: bool = False,
        export_only: bool = False,
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
        pti = parse_transcript_id
        if pti is None:
            pti = TranscriptID.parser_from_re(
                config.transcript_id_regex,
            )
        exporter = export_orthologs.OrthologExporter(
            sim,
            pti,
            False,
            debug=True,
            consistent_strands=True,
            allow_inconsistent=True,
            jobs=config.jobs,
        )
        component_paths = exporter.by_component(export_dir, order="after")
        if not export_only:
            all_ideal_path = export_dir / "all_ideal.fasta"
            with open(all_ideal_path, "w") as all_ideal:
                for path in component_paths.values():
                    with open(path, "r") as component_fasta:
                        all_ideal.write(component_fasta.read())
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
                    node_to_ccc=exporter.node_to_component_component
                )
                # if stats is None:
                #     from IPython import embed; embed()
                with open(search_dir / "stats", "w") as stats_file:
                    json.dump(stats._asdict(), stats_file)

def main():
    _, args = build_parser().get_arguments()
    configs = [config_module.RNACliqueConfig.yaml_load(c) for c in args.configs]
    args.export_output_dir.mkdir(exist_ok=True)
    parse_transcript_id = None
    if args.transcript_id_regex is not None:
        parse_transcript_id = TranscriptID.parser_from_re(
            args.transcript_id_regex
        )
    try:
        export_and_search(
            configs,
            args.export_output_dir,
            args.queries,
            parse_transcript_id,
            args.jobs,
            args.resolve_name_conflicts,
            args.export_only
        )
    except NameConflictError as e:
        eprint(e)
        eprint(("Cannot continue. Provide the -r option to try automatic "
                "resolution."))
        sys.exit(1)

if __name__ == "__main__":
    main()
