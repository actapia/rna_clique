import argparse
import json
import search_ideal_components
import export_orthologs
from make_subset import multi_glob
from filtered_distance import SampleSimilarity

from collections import defaultdict
from tqdm import tqdm

from pathlib import Path

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-A", "--analyses", type=Path, nargs="+", required=True)
    parser.add_argument("-O", "--out-root", type=Path, required=True)
    parser.add_argument("-j", "--jobs", type=int, default=1)
    parser.add_argument("-X", "--export-only", action="store_true")
    parser.add_argument("-Q", "--queries", nargs="+", type=Path, required=True)
    return parser.parse_args()

def select_out_dir_names(analyses):
    out_to_analysis = {}
    to_fix = defaultdict(set)
    for a in analyses:
        try:
            new = {(a,)*2, out_to_analysis[a.name]}
            to_fix[a.name] |= new
        except KeyError:
            out_to_analysis[a.name] = (a, a)
    print(out_to_analysis, to_fix)
    while to_fix:
        new_to_fix = {}
        for k, s in to_fix.items():
            del out_to_analysis[k]
            for analysis, ancestor in s:
                if not ancestor.parent.name:
                    raise ValueError("Could not resolve output name conflicts.")
                try:
                    new = {
                        (analysis, ancestor.parent),
                        out_to_analysis[ancestor.parent.name]
                    }
                    new_to_fix[ancestor.parent.name] |= new
                except KeyError:
                    out_to_analysis[ancestor.parent.name] = (
                        analysis,
                        ancestor.parent
                    )
        to_fix = new_to_fix
    return {v[0]: k for (k, v) in out_to_analysis.items()}

def main():
    args = handle_arguments()

    out_dir_names = select_out_dir_names(args.analyses)
    for analysis in tqdm(args.analyses):
        out_dir = args.out_root / out_dir_names[analysis]
        out_dir.mkdir(exist_ok=True)
        export_dir = out_dir / "export"
        export_dir.mkdir(exist_ok=True)
        graph_path =  analysis / "graph.pkl"
        comparison_paths = list(
            multi_glob(analysis / "od2", ["*.pkl", "*.h5"])
        )
        sim = SampleSimilarity.from_filenames(
            graph_path,
            tqdm(comparison_paths),
            store_dfs=True
        )
        exporter = export_orthologs.OrthologExporter(
            sim,
            export_orthologs.default_gene_re,
            False,
            debug=True,
            consistent_strands=True,
            allow_inconsistent=True,
            jobs=args.jobs,
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
                    graph_loc=graph_path,
                    comparisons_loc=comparison_paths,
                    exported=all_ideal_path,
                    db_cache_loc=db_cache,
                    out_dir=search_dir,
                    query=query,
                    debug=True,
                    sample_count=None,
                    clean=True,
                    merge_sams=True,
                    jobs=args.jobs,
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
