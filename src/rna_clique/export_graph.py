import json
import io
import sys
import pickle

import networkx as nx


from pathlib import Path
from contextlib import ExitStack

from . import config as config_module
from . import app
from .app import set_except_hook, eprint, get_format_from_extension

def write_cytoscape(graph: nx.Graph, out_file: io.TextIOBase):
    """Export the given graph as a Cytoscape.js JSON file.

    Parameters:
        graph:    The graph to export.
        out_file: Output file-like object.
    """
    json.dump(nx.cytoscape_data(graph), out_file)


def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager(
        description="Export a gene matches graph to another format."
    )
    arg_config.expose_fields_with_default_aliases(
        "graph",
        required=True
    )
    arg_config.expose_config_field(
        "output_dir",
        aliases=["--analysis-root", "--rna-clique-output-dir", "-A"]        
    )
    arg_config.add_argument(
        "--export-out",
        "-x",
        type=Path,
        help="Path to which to export the graph."
    )
    arg_config.add_argument(
        "-f",
        "--format",
        choices=writers,
        default={
            ("export_out",): get_format_from_extension(extension_to_format)
        },
        help="Format for writing graph.",
        required=True
    )    
    return arg_config

extension_to_format = {
    "cyjs": "cytoscape",
    "graphml": "graphml",
    "dot": "graphviz",
    "gv": "graphviz"
}

# type_name = {
#     "cytoscape": "Cytoscape.js JSON",
#     "graphml": "GraphML",
#     "graphviz": "Graphviz"
# }

def text_io_wrapped(fun):
    def inner(graph, f):
        return fun(graph, io.TextIOWrapper(f))
    return inner

writers = {
    "cytoscape": text_io_wrapped(write_cytoscape),
    "graphml": nx.write_graphml,
    "graphviz": text_io_wrapped(nx.drawing.nx_pydot.write_dot)
}

def main():
    with set_except_hook():
        try:
            _, args, config = build_parser().get_arguments_and_config()
        except app.UnrecognizedFileExtensionError:
            eprint(
                "Could not determine file format automatically from "
                "extension for output graph filename. Please change the "
                "file extension or specify a format explicitly with the "
                "--format/-f option.\n"
            )
            raise
    with set_except_hook(args.verbose):
        with open(config.graph, "rb") as graph_pickle:
            graph = pickle.load(graph_pickle)
        with ExitStack() as stack:
            if args.export_out:
                f = open(args.export_out, "wb")
                stack.push(f)
            else:
                f = sys.stdout.buffer
            writers[args.format](graph, f)

if __name__ == "__main__":
    main()
