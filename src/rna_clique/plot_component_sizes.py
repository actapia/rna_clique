import pickle

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path
from typing import Optional

from . import config as config_module
from .graph import component_subgraphs
from .gene_matches_tables import read_table, get_table_files
from .app import set_except_hook, eprint

def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager(
        description="Produce visualizations based on a gene matches graph.",
    )
    arg_config.expose_fields_with_default_aliases(
        "graph",
        required=True
    )
    arg_config.expose_fields_with_default_aliases(
        "tables_dir"
    )
    arg_config.expose_config_field(
        "output_dir",
        aliases=["--analysis-root", "--rna-clique-output-dir", "-A"]        
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
    # arg_config.add_argument(
    #     "-G",
    #     "--graphviz",
    #     type=Path,
    #     help="output path for graphviz (dot) representation"
    # )
    # arg_config.add_argument(
    #     "-x",
    #     "--export",
    #     nargs="+",
    #     type=Path,
    #     help="output paths for export to {}.".format(
    #         " or ".join(type_name.values())
    #     )
    # )
    arg_config.add_argument(
        "--statistics",
        nargs="?",
        choices=["h", "m"],
        const="h",
        help=("print statistics in the desired format (human or "
              "machine-readable)")
    )
    return arg_config

stat_labels = [
    "Samples",
    "Total components",
    "Components >= samples",
    "Ideal components"
]

def component_hist(data: list[int], highlight: Optional[int] = None):
    """Draw a histogram for the given data.

    This function is made for drawing histograms for data in small integer
    ranges, so the bin size is always 1.

    One bar of the histogram can be highlighted in a different color. This could
    represent, for example, the bar corresponding to the number of samples
    present in the analysis.

    Parameters:
        data (list):     Data from which to compute frequency distribution.
        highlight (int): Value frequency to highlight in a different color.
    """
    hist, bins = np.histogram(
        data,
        bins=range(1, max(data) + 2)
    )
    # embed()
    plt.bar(bins[:-1], hist, width=np.diff(bins), color="C0", align="center")
    plt.bar(
        bins[highlight - 1],
        hist[highlight - 1],
        width=np.diff(bins),
        color="C1",
        align="center"
    )

class SampleCountError(ValueError):
    pass

def count_samples(config: config_module.RNACliqueConfig) -> int:
    """Determine the number of samples used in an analysis.

    This function prefers faster ways of determining the number of samples if
    they are available. If the number of samples could not be determined from
    the command-line arguments and configuration, then a SampleCountError will
    be raised.

    Parameters:
        config: RNACliqueConfig for the analysis.

    Returns:
        The number of samples used in the analysis.

    """
    if config.path_to_sample is not None:
        return len(config.path_to_sample)
    if config.input_dirs is not None:
        return len(config.input_dirs)
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
    raise SampleCountError(
        ("Could not determine number of samples from arguments "
         "and configuration.")
    )

def main():
    with set_except_hook():
        _, args, config = build_parser().get_arguments_and_config()
    with set_except_hook(config.verbose):
        try:
            samples = count_samples(config)
        except SampleCountError:
            eprint("Count not determine number of samples for the analysis. "
                   "Please provide a config file with the path_to_sample "
                   "attribute set, or provide the tables_dir setting.\n")
            raise
        with open(config.graph, "rb") as f:
            graph = pickle.load(f)
        components = list(component_subgraphs(graph))
        # embed()
        component_sizes = [len(c) for c in components]
        sample_counts = [len(set(t[0] for t in c)) for c in components]
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
            density = [
                2*len(c.edges)/(l*(l-1))
                for (c, l) in zip(components, component_sizes)
            ]        
            eprint("Making density plot.")
            sns.set_style("whitegrid")
            kde = sns.kdeplot(np.array(sorted(density)))
            fig = kde.get_figure()
            fig.savefig(args.density_plot)
        # if args.graphviz:
        #     eprint("Making graphviz file.")
        #     nx.drawing.nx_pydot.write_dot(graph, args.graphviz)
        # if args.export is not None:
        #     for exp in args.export:
        #         type_ = extension_to_type[exp.suffix]
        #         eprint("Exporting to {}.".format(type_name[type_]))
        #         exporters[type_](graph, exp)
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
