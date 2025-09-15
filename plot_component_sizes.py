import argparse
import pickle
import json
import config as config_module

from pathlib import Path

import numpy as np
import networkx as nx

import matplotlib.pyplot as plt
import seaborn as sns

from find_homologs import eprint
from graph import component_subgraphs
from gene_matches_tables import read_table, get_table_files

def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager()
    arg_config.expose_fields_with_default_aliases(
        "graph",
        required=True
    )
    arg_config.expose_fields_with_default_aliases(
        "top_genes_dir",
        "tables_dir"
    )
    arg_config.expose_config_field(
        "output_dir",
        aliases=["--analysis-root", "--rna-clique-output-dir", "-A"]        
    )
    arg_config.add_argument(
        "--samples",
        type=int,
        help="number of samples"
    )
    arg_config.add_argument(
        "-s",
        "--size-plot",
        type=Path,
        help="output path for component size histogram"
    )
    arg_config.add_argument(
        "-S",
        "--sample-plot",
        type=Path,
        help="output path for represented sample count plot"
    )
    arg_config.add_argument(
        "-r",
        "--ratio-plot",
        type=Path,
        help="output path for KDE of represented sample count / component size"
    )
    arg_config.add_argument(
        "-d",
        "--density-plot",
        type=Path,
        help="output path for KDE of component density"
    )
    arg_config.add_argument(
        "-G",
        "--graphviz",
        type=Path,
        help="output path for graphviz (dot) representation"
    )
    arg_config.add_argument(
        "-x",
        "--export",
        nargs="+",
        type=Path,
        help="output paths for export to {}.".format(
            " or ".join(type_name.values())
        )
    )
    arg_config.add_argument(
        "--statistics",
        nargs="?",
        choices=["h", "m"],
        const="h",
        help=("print statistics in the desired format (human or "
              "machine-readable)")
    )
    return arg_config

def export_cytoscape(graph : nx.Graph, out_file : Path):
    """Export the given graph as a Cytoscape JSON file."""
    with open(out_file, "w") as jf:
        json.dump(nx.cytoscape_data(graph), jf)

extension_to_type = {
    ".cyjs": "cytoscape",
    ".graphml": "graphml"
}

type_name = {
    "cytoscape": "Cytoscape JSON",
    "graphml": "GraphML"
}

exporters = {
    "cytoscape": export_cytoscape,
    "graphml": nx.write_graphml
}

stat_labels = [
    "Samples",
    "Total components",
    "Components >= samples",
    "Ideal components"
]

def component_hist(data: list[int], samples: int):
    hist, bins = np.histogram(
        data,
        bins=range(1, max(data) + 2)
    )
    # embed()
    plt.bar(bins[:-1], hist, width=np.diff(bins), color="C0", align="center")
    plt.bar(
        bins[samples - 1],
        hist[samples - 1],
        width=np.diff(bins),
        color="C1",
        align="center"
    )

def count_samples(args, config):
    if args.samples is not None:
        return args.samples
    if config.path_to_sample is not None:
        return len(config.path_to_sample)
    if config.input_dirs is not None:
        return len(config.input_dirs)
    if config.top_genes_dir is not None:
        return sum(1 for _ in config.top_genes_dir.glob("*.fasta"))
    if config.tables_dir is not None:
        return len(
            set.union(
                *(
                    {df.iloc[0][f"{s}sample"] for s in "sq"} for df in
                    map(
                        read_table,
                        get_table_files(config.tables_dir)
                    )
                )
            )
        )
    raise ValueError(
        ("Could not determine number of samples from arguments "
         "and configuration.")
    )
                        

def main():
    _, args, config = build_parser().get_arguments_and_config()
    samples = count_samples(args, config)
    with open(config.graph, "rb") as f:
        graph = pickle.load(f)
    components = list(component_subgraphs(graph))
    # embed()
    component_sizes = [len(c) for c in components]
    sample_counts = [len(set(t[0] for t in c)) for c in components]
    density = [
        2*len(c.edges)/(l*(l-1)) for (c, l) in zip(components, component_sizes)
    ]
    if args.size_plot:
        eprint("Making size plot.")
        # cs_counter = Counter(component_sizes)
        # max_size = max(cs_counter)
        # size_counts = [cs_counter[k] for k in range(1, max_size)]
        component_hist(component_sizes, samples)
        #plt.hist(component_sizes, bins=range(1, max_size))
        plt.xlabel("Component size")
        plt.ylabel("Frequency")
        plt.savefig(args.size_plot)
    plt.clf()
    if args.sample_plot:
        eprint("Making sample plot.")
        # sc_counter = Counter(sample_counts)
        # max_size = max(sc_counter)
        # size_counts = [sc_counter[k] for k in range(1, max_size)]
        component_hist(sample_counts, samples)
        plt.xlabel("Sample count")
        plt.ylabel("Frequency")
        plt.savefig(args.sample_plot)
    plt.clf()
    if args.ratio_plot:
        eprint("Making ratio plot.")
        ratios = [a/b for (a, b) in zip(sample_counts, component_sizes)]
        sns.set_style("whitegrid")
        kde = sns.kdeplot(np.array(sorted(ratios)))
        fig = kde.get_figure()
        fig.savefig(args.ratio_plot)
    plt.clf()
    if args.density_plot:
        eprint("Making density plot.")
        sns.set_style("whitegrid")
        kde = sns.kdeplot(np.array(sorted(density)))
        fig = kde.get_figure()
        fig.savefig(args.density_plot)
    if args.graphviz:
        eprint("Making graphviz file.")
        nx.drawing.nx_pydot.write_dot(graph, args.graphviz)
    if args.export is not None:
        for exp in args.export:
            type_ = extension_to_type[exp.suffix]
            eprint("Exporting to {}.".format(type_name[type_]))
            exporters[type_](graph, exp)
    if args.statistics:
        eprint("Computing statistics.")
        ideal = sum(
            1 for (c, s, g) in zip(
                components,
                component_sizes,
                sample_counts
            )
            if len(c) == samples and g == samples and
            2*len(c.edges) == s*(s-1)
        )
        gt_samples = sum(1 for c in components if len(c) >= samples)
        stats = [samples, len(components), gt_samples, ideal]
        if args.statistics == "h":
            for label, stat in zip(stat_labels, stats):
                print(label + ":", stat)
        else:
            print(*stats)
        
        #embed()
        
if __name__ == "__main__":
    main()
