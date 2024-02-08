import argparse
import pickle
import json

from fractions import Fraction

from collections import Counter
from pathlib import Path
from typing import Iterator

import numpy as np
import networkx as nx

import matplotlib.pyplot as plt
import seaborn as sns

from find_homologs import eprint

from IPython import embed

def handle_arguments():
    parser = argparse.ArgumentParser(
        description=("generate exports, visualizations, or statistics for a "
                     "gene matches graph")
    )
    #parser.add_argument("-i", "--inputs", nargs="+", type=Path, required=True)
    parser.add_argument(
        "graph",
        type=Path,
        help="path to the gene matches graph pickle"
    )
    parser.add_argument(
        "-s",
        "--size-plot",
        type=Path,
        help="output path for component size histogram"
    )
    parser.add_argument(
        "-S",
        "--sample-plot",
        type=Path,
        help="output path for represented sample count plot"
    )
    parser.add_argument(
        "-r",
        "--ratio-plot",
        type=Path,
        help="output path for KDE of represented sample count / component size"
    )
    parser.add_argument(
        "-d",
        "--density-plot",
        type=Path,
        help="output path for KDE of component density"
    )
    parser.add_argument(
        "-g",
        "--graphviz",
        type=Path,
        help="output path for graphviz (dot) representation"
    )
    parser.add_argument(
        "-x",
        "--export",
        nargs="+",
        type=Path,
        help="output paths for export to {}.".format(
            " or ".join(type_name.values())
        )
    )
    parser.add_argument(
        "--samples",
        type=int,
        help="number of samples in the analysis"
    )
    parser.add_argument(
        "--statistics",
        nargs="?",
        choices=["h", "m"],
        const="h",
        help=("print statistics in the desired format (human or "
              "machine-readable)")
    )
    return parser.parse_args()

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


def component_subgraphs(g : nx.Graph) -> Iterator[nx.Graph]:
    """Yields the connected components of the given graph as subgraphs."""
    for c in nx.connected_components(g):
        yield g.subgraph(c)

def main():
    args = handle_arguments()
    with open(args.graph, "rb") as f:
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
        hist, bins = np.histogram(
            component_sizes,
            bins=range(1, max(component_sizes) + 2)
        )
        #embed()
        plt.bar(bins[:-1], hist, width=np.diff(bins), color="C0", align="center")
        plt.bar(
            bins[args.samples-1],
            hist[args.samples-1],
            width=np.diff(bins),
            color="C1",
            align="center"
        )
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
        hist, bins = np.histogram(
            sample_counts,
            bins=range(1, max(sample_counts) + 2)
        )
        # embed()
        plt.bar(bins[:-1], hist, width=np.diff(bins), color="C0", align="center")
        plt.bar(
            bins[args.samples-1],
            hist[args.samples-1],
            width=np.diff(bins),
            color="C1",
            align="center"
        )
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
        if args.samples is None:
            args.samples = len(set(t[0] for t in graph.nodes))
        ideal = sum(
            1 for (c, s, g) in zip(
                components,
                component_sizes,
                sample_counts
            )
                    if len(c) == args.samples and g == args.samples and \
            2*len(c.edges) == s*(s-1)
        )
        gt_samples = sum(1 for c in components if len(c) >= args.samples)
        stats = [args.samples, len(components), gt_samples, ideal]
        if args.statistics == "h":
            for label, stat in zip(stat_labels, stats):
                print(label + ":", stat)
        else:
            print(*stats)
        
        #embed()
        
if __name__ == "__main__":
    main()
