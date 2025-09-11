import argparse
import json
import sys
import re
import search_ideal_components
import export_orthologs
import config as config_module
from filtered_distance import SampleSimilarity
from name_conflict_resolver import NameConflictResolver
from find_homologs import eprint
from gene_matches_tables import get_table_files

from collections import defaultdict, Counter
from tqdm import tqdm

from transcripts import default_gene_re, TranscriptID

from pathlib import Path

def handle_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument("-A", "--analyses", type=Path, nargs="+", required=True)
    parser.add_argument("-C", "--configs", type=Path, nargs="+", required=True)
    parser.add_argument("-r", "--resolve-name-conflicts", action="store_true")
    parser.add_argument("-X", "--export-output-dir", type=Path, required=True)
    parser.add_argument("-j", "--jobs", type=int)
    parser.add_argument("-x", "--export-only", action="store_true")
    parser.add_argument("-Q", "--queries", nargs="+", type=Path, required=True)
    parser.add_argument(
        "--transcript-id-regex",
        "-p",
        "--pattern",
        type=re.compile,
        help="Python regex for parsing sequence IDs",
        default=default_gene_re
    )
    return parser.parse_args()

def get_analysis_name(config):
    name = config.name
    if not name:
        try:
            out_dir = config.output_dir
        except AttributeError:
            raise ValueError("Could not determine analysis name.")
        name = out_dir.name
    return name

# def select_out_dir_names(analyses, resolve=True):
#     out_to_analysis = {}
#     to_fix = defaultdict(set)
#     for a in analyses:
#         try:
#             new = {(a,)*2, out_to_analysis[get_analysis_name(a)]}
#             to_fix[a.name] |= new
#         except KeyError:
#             out_to_analysis[a.name] = (a, a)
#     if to_fix and not resolve:
#         raise ValueError("Could not resolve output name conflicts.")
#     while to_fix:
#         new_to_fix = {}
#         for k, s in to_fix.items():
#             del out_to_analysis[k]
#             for analysis, ancestor in s:
#                 if not ancestor.parent.name:
#                     raise ValueError("Could not resolve output name conflicts.")
#                 try:
#                     new = {
#                         (analysis, ancestor.parent),
#                         out_to_analysis[ancestor.parent.name]
#                     }
#                     new_to_fix[ancestor.parent.name] |= new
#                 except KeyError:
#                     out_to_analysis[ancestor.parent.name] = (
#                         analysis,
#                         ancestor.parent
#                     )
#         to_fix = new_to_fix
#     return {v[0]: k for (k, v) in out_to_analysis.items()}

def main():
    args = handle_arguments()
    configs = [config_module.RNACliqueConfig.yaml_load(c) for c in args.configs]
    out_names = [args.export_output_dir / get_analysis_name(c) for c in configs]
    counts = Counter(out_names)
    try:
        multiple = next(x for (x, c) in counts.items() if c > 1).name
        eprint(f"Multiple analyses named {multiple!r}.")
        if not args.resolve_name_conflicts:
             eprint(("Cannot continue. Provide the -r option to try automatic "
                     "resolution."))
             sys.exit(1)
        else:
            eprint("Resolving conflicts automatically.")
    except StopIteration:
        pass
    resolver = NameConflictResolver.from_keys(
        range(len(out_names)),
        out_names.__getitem__
    )
    for config, (_, (_, out_dir)) in zip(configs, resolver.resolve()):
        if config.transcript_id_regex is None:
            config.transcript_id_regex = args.transcript_id_regex
        if args.jobs is not None:
            config.jobs = args.jobs
        out_dir.mkdir(exist_ok=True)
        export_dir = out_dir / "export"
        export_dir.mkdir(exist_ok=True)
        comparison_paths = list(get_table_files(config.tables_dir))
        sim = SampleSimilarity.from_filenames(
            config.graph,
            tqdm(comparison_paths),
            store_dfs=True
        )
        exporter = export_orthologs.OrthologExporter(
            sim,
            config.transcript_id_regex,
            False,
            debug=True,
            consistent_strands=True,
            allow_inconsistent=True,
            jobs=config.jobs,
        )
        component_paths = exporter.by_component(export_dir, order="after")
        if not args.export_only:
            all_ideal_path = export_dir / "all_ideal.fasta"
            with open(all_ideal_path, "w") as all_ideal:
                for path in component_paths.values():
                    with open(path, "r") as component_fasta:
                        all_ideal.write(component_fasta.read())
            db_cache = export_dir / "db_cache"
            db_cache.mkdir(exist_ok=True)
            for query in args.queries:
                search_dir = out_dir / ("search_" + query.stem)
                search_dir.mkdir(exist_ok=True)
                stats = search_ideal_components.search(
                    graph_loc=config.graph,
                    comparisons_loc=comparison_paths,
                    exported=all_ideal_path,
                    db_cache_loc=db_cache,
                    out_dir=search_dir,
                    query=query,
                    debug=True,
                    sample_count=None,
                    clean=True,
                    merge_sams=True,
                    parse_transcript_ids=TranscriptID.parse_from_re(
                        config.transcript_id_regex,
                    ),
                    jobs=config.jobs,
                    sim=sim,
                    strand_graph_out=(
                        exporter.strand_graph,
                        exporter.node_to_component_component
                    )
                )
                if stats is None:
                    from IPython import embed; embed()
                with open(search_dir / "stats", "w") as stats_file:
                    json.dump(stats._asdict(), stats_file)

if __name__ == "__main__":
    main()
