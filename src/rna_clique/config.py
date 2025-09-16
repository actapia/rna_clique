import importlib
import multiprocessing
import argparse
import dataclasses
import sys
import typing
import json
import functools
import math
import itertools
import re
import pprint
import copy

import networkx as nx

from datetime import datetime
from typing import Optional
from collections import deque, ChainMap
from collections.abc import Collection

from .transcripts import default_gene_re
from .marshalling_dataclass import marshalling_dataclass, marshalling_field

from pathlib import Path
from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

def get_version():
    version = None
    try:
        import git
        try:
            repo = git.Repo(
                Path(__file__).parent,
                search_parent_directories=True
            )
            version = repo.head.object.hexsha
            changed= [i.a_path for i in repo.index.diff(None)]
            if changed:
                version += "+dev"
        except git.InvalidGitRepositoryError:
            pass
    except ImportError:
        pass
    if version is None:
        try:
            rnac = importlib.import_module(__package__)
            files = {
                str(Path(rnac.__file__).parent / f)
                for f in importlib.metadata.files(rnac.__name__)
            }
            if __file__ in files:
                version = importlib.metadata.version(__package__)
        except ValueError:
            pass
    if version is None:
        version = "dev"
    return version

@marshalling_dataclass(optional=True)
class RNACliqueConfig:
    config_version: str = marshalling_field(default="0.0.1", metadata={
        "description": "Version of the configuration schema used."})
    title: Optional[str] = marshalling_field(metadata={
        "description": "Name to assign to the analysis."
    })
    input_dirs: Optional[list[Path]] = marshalling_field(list[str], metadata={
        "description": "Directories containing the transcript FASTA files."})
    top_genes_dir: Optional[Path] = marshalling_field(str, metadata={
        "description": "Output directory containing top n genes by coverage."})
    tables_dir: Optional[Path] = marshalling_field(str, metadata={
        "description": "Output directory containing gene matches tables."})
    cache_dir: Optional[Path] = marshalling_field(str, metadata={
        "description": "Intermediate directory containing BLAST DB caches."})
    output_dir: Optional[Path] = marshalling_field(str, metadata={
        "description": "Output directory root."})
    graph: Optional[Path] = marshalling_field(str, metadata={
        "description": "Output gene matches graph."})
    top_genes: Optional[int] = marshalling_field(metadata={
        "description": "Number of top genes by k-mer coverate to select."})
    transcripts_name: Optional[str] = marshalling_field(
        default="transcripts.fasta",
        metadata={
            "description": "Name of transcripts files in input directories."
        }
    )
    top_matches: Optional[int] = marshalling_field(default=1, metadata={
        "description": "Threshold for counting a match between two genes."})
    evalue: Optional[float] = marshalling_field(default=1e-99, metadata={
        "description": "e-value threshold to use for BLASTn searches."})
    keep_all: Optional[bool] = marshalling_field(default=True, metadata={
        "description": "Keep all matches between genes in the case of ties."})
    jobs: Optional[int] = marshalling_field(
        default=multiprocessing.cpu_count() - 1,
        metadata={
            "description": "Number of parallel jobs to use."
        }
    )
    transcript_id_regex: Optional[re.Pattern] = marshalling_field(
        lambda x: x.pattern,
        unmarshal=re.compile,
        default=default_gene_re,
        metadata={
            "description": "Python regex to use for parsing transcript IDs."
        }
    )
    path_to_sample: Optional[dict[Path, str]] = marshalling_field(
        dict[str, str],
        metadata={
            "description": "Mapping from paths to sample names."
        }
    )
    matrix: Optional[Path] = marshalling_field(str, metadata={
        "description": "Output distance matrix location."})
    finished: Optional[datetime] = marshalling_field(
        marshal=datetime.isoformat,
        unmarshal=datetime.fromisoformat,
        metadata={
            "description": ("When the last analysis associated with this "
                            "config file finished.")
        }
    )
    version: str = marshalling_field(default=get_version(), metadata={
        "description": "Version of RNA-clique used to create this analysis."})
    subset_of: Optional[Path] = marshalling_field(
        str,
        metadata={
            "description": "Path to analysis of which this is a subset."
        }
    )

    @classmethod
    def yaml_load(cls, path: Path):
        with open(path, "r") as f:
            return cls.from_marshalled_representation(
                load(
                    f,
                    Loader=Loader
                )
            )

    def mark_finish(self):
        self.finished = datetime.now()

    def yaml_save(self, path: Path):
        with open(path, "w") as f:
            f.write(dump(self.marshal(hide_none=True), Dumper=Dumper))

    @classmethod
    def from_args(cls, args):
        if args.input_config is not None:
            config = cls.yaml_load(args.input_config)
        else:
            config = cls()
        for field in cls.__dataclass_fields__:
            arg_value = None
            try:
                arg_value = getattr(args, field)
            except AttributeError:
                pass
            if getattr(config, field) is None or arg_value:
                setattr(config, field, arg_value)
        return config

# def get_args_and_config(parser, ins=None, outs=None):
#     if ins is None:
#         ins = {}
#     if outs is None:
#         outs = {}
#     args = parser.parse_args()
#     process_out_dir_files(args, ins, outs)
#     config = RNACliqueConfig.from_args(args)
#     process_config(parser, config, ins, outs)
#     return args, config

# def process_out_dir_files(x, ins, outs):
#     for arg, f in outs.items():
#         if getattr(x, arg.dest, None) is None:
#             setattr(x, arg.dest, x.output_dir / f)
#     for arg, f in ins.items():
#         in_file = x.output_dir / f
#         if in_file.exists() and getattr/(x, arg.dest, None) is None:
#             setattr(x, arg.dest, x.output_dir / f)

# def detect_missing(parser, config, args, make_message):
#     try:
#         missing, _ = next(
#             x for x in args if x.required and getattr(config, x.dest) is None
#         )
#         parser.error(make_message(missing))
#     except StopIteration:
#         pass    
            
# def process_config(
#         parser,
#         config,
#         ins,
#         outs,
# ):
#     process_out_dir_files(config, ins, outs)
#     detect_missing(
#         parser,
#         config,
#         outs,
#         lambda x: "Must provide --output-dir or --{}".format(
#             x.replace("_", "_")
#         )
#     )
#     detect_missing(
#         parser,
#         config,
#         ins,
#         lambda x: ("Must provide --{}, or {} must be present in "
#                    "--output-dir").format(
#                        x.replace("_","-"),
#                        ins[x],
#                    )
#     )

# def typing_to_arg_collection_and_element_type(t):
#     fields = typing.get_args(t)
#     type_ = typing.get_origin(t)
#     if type_ is not None:
#         if isinstance(type_, Collection):
#             if len(fields) == 1:
#                 pass

def detect_nargs(t, multi_type="+"):
    fields = typing.get_args(t)
    type_ = typing.get_origin(t)
    if type_ is not None:
        if isinstance(type_, type):
            if issubclass(type_, tuple) or issubclass(type_, Collection):
                return multi_type
        elif type_ is typing.Union:
            non_none = [f for f in fields if f is not type(None)]    
            if len(fields) == 2 and len(non_none) == 1:
                # Handle "Optional" types.
                non_none = non_none[0]
                return detect_nargs(non_none, multi_type=multi_type)
    return None

def get_edge_priority(graph, edge):
    return graph.edges[edge]["priority"]

def update_node(graph, origin):
    changed = False
    priority_key = functools.partial(get_edge_priority, graph)
    for priority, edge_group in itertools.groupby(
            sorted(
                graph.in_edges(origin),
                key=priority_key
            ),
            priority_key
        ):
            if priority > graph.nodes[origin]["priority"]:
                break
            edge = next(iter(edge_group))
            kwargs = {
                x: graph.nodes[x]["value"]
                for x in graph.edges[edge]["args"]
            }
            if all(kwargs.values()):
                new_value = graph.edges[edge]["function"](**kwargs)
                if new_value is not None:
                    if new_value != graph.nodes[origin]["value"]:
                        graph.nodes[origin]["value"] = new_value
                        changed = True
                    graph.nodes[origin]["priority"] = graph.edges[
                        edge
                    ]["priority"]
                    break
    return changed
    

def propagate_defaults(graph):
    con = nx.condensation(graph)
    # print(
    #     "condensation has {} nodes, {} edges".format(
    #         len(list(con.nodes)),
    #         len(list(con.edges))
    #     )
    # )
    
    for scc in nx.topological_sort(con):
        subgraph = graph.subgraph(con.nodes[scc]["members"])
        for pred in con.predecessors(scc):
            # print("pred is", con.nodes[pred]["members"])
            # print("scc is", con.nodes[scc]["members"])
            for node in nx.node_boundary(
                    graph,
                    con.nodes[pred]["members"],
                    con.nodes[scc]["members"]
            ):
                
                update_node(graph, node)
        st = deque(
            (n, True)
            for n in subgraph.nodes if subgraph.nodes[n]["value"] is not None
        )
        while st:
            origin, changed = st.pop()
            changed |= update_node(subgraph, origin)
            if changed:
                for edge in subgraph.out_edges(origin):
                    dest = edge[1]
                    if subgraph.edges[edge]["priority"] <= \
                       subgraph.nodes[dest]["priority"]:
                        st.append(dest)

class DefaultDiGraph(nx.DiGraph):
    def __init__(self, *args, node_defaults={}, edge_defaults={}, **kwargs):
        super().__init__(*args, **kwargs)
        self.node_defaults = node_defaults
        self.edge_defaults = edge_defaults

    def node_attr_dict_factory(self):
        return dict(self.node_defaults)

    def edge_attr_dict_factory(self):
        return dict(self.edge_defaults)

class MultiArgConstAction(argparse.Action):
    def _get_values(self, values):
        if not values:
            return self.const
        else:
            return values
        
    def __call__(self, parser, namespace, values, option_string):
        setattr(namespace, self.dest, self._get_values(values))

class UnmarshallingStoreAction(MultiArgConstAction):
    def __init__(self, *args, unmarshal=lambda x: x, **kwargs):        
        super().__init__(*args, **kwargs)
        self.unmarshal = unmarshal

    def _get_values(self, values):
        return self.unmarshal(super()._get_values(values))


class ConfigArgumentManager:
    default_node = object()
    _type_to_action = {
        bool: "store_true",
        Optional[bool]: "store_true"
    }
    
    def __init__(
            self,
            config_class,
            *args,
            default_deps={},
            accept_input=True,
            show_config=True,
            show_config_format=True,
            #use_post_type=False,
            type_to_action=None,
            multi_nargs_type="+",
            **kwargs
    ):
        if type_to_action is None:
            self._type_to_action = dict(type(self)._type_to_action)
        self._config_class = config_class
        self.parser = argparse.ArgumentParser(*args, **kwargs)
        self._default_graph = DefaultDiGraph(
            node_defaults={
                "value": None,
                "priority": math.inf,
                "required": False,
            }
        )
        self._option_string_to_dest = {}
        self._dest_to_option_strings = {}#defaultdict(list)
        for name, f in self._config_class.__dataclass_fields__.items():
            self._default_graph.add_node(name, value=f.default)
        for dest, deps in default_deps.items():
            self.add_defaults(dest, deps)
        if accept_input:
            kw = {}
            if accept_input is not True:
                kw = {"aliases": accept_input}
            self.add_config_argument(**kw)
        if show_config:
            kw = {}
            if show_config is not True:
                kw = {"aliases": show_config}
            self.add_show_config_argument(**kw)
        if show_config_format:
            kw = {}
            if show_config_format is not True:
                kw = {"aliases": show_config_format}
            self.add_show_config_format_argument(**kw)
        self._multi_nargs_type = multi_nargs_type
        # self._post_type = defaultdict(lambda: lambda x: x)
        # self._use_post_type = use_post_type

    def add_config_argument(self, aliases=["--input-config", "-c", "-c1"]):
        self.add_argument(
            *aliases,
            dest="input_config",
            type=Path
        )

    def add_show_config_argument(self, aliases=["--show-config"]):
        self.add_argument(
            *aliases,
            dest="_config_show_config",
            nargs="*",
            action=MultiArgConstAction,
            const=["config"],
            choices=["original_args", "args", "config"]
        )

    @classmethod
    def select_default_config_format(cls, _config_show_config):
        if _config_show_config == ["config"]:
            return "yaml"
        else:
            return "dict"

    def add_show_config_format_argument(self, aliases=["--show-config-format"]):
        self.add_argument(
            *aliases,
            dest="_config_show_config_format",
            choices=["dict", "yaml", "json"],
            default={
                (
                    "_config_show_config",
                ): type(self).select_default_config_format
            }
        )

    def add_defaults(self, option, default, next_priority=None):
        #print("\tAdd defaults",option,default)
        dest = self._option_string_to_dest.get(option, option)
        if next_priority is None:
            next_priority = max(
                (
                    self._default_graph.edges[e]["priority"]
                    for e in self._default_graph.in_edges(dest)
                ),
                default=-1
            ) + 1
        try:
            for priority, (parents, fun) in enumerate(
                    default.items(),
                    next_priority
            ):
                if parents:
                    #print("\t\tAdding for",parents)
                    for parent in parents:
                        self._default_graph.add_edge(
                            parent, dest,
                            priority=priority,
                            function=fun,
                            args=parents
                        )                        
                else:
                    self._default_graph.add_edge(
                        self.default_node, dest,
                        priority=priority,
                        function=lambda: fun,
                        args=(),
                    )
        except AttributeError:
            self._default_graph.add_edge(
                self.default_node, dest,
                priority=next_priority,
                function=lambda: default,
                args=(),
            )

    def set_defaults(self, option, default):
        dest = self._option_string_to_dest.get(option, option)
        self._default_graph.remove_edges_from(
            list(self._default_graph.in_edges(dest))
        )
        self.add_defaults(option, default)

    def add_argument(
            self,
            *args,
            required=False,
            default=None,
            #post_type=None,
            **kwargs
    ):
        try:
            action = self.parser.add_argument(*args, **kwargs)
        except ValueError:
            from IPython import embed; embed()
        self._default_graph.add_node(action.dest, required=required)
        self.add_defaults(action.dest, default)
        self._option_string_to_dest |= {
            s: action.dest for s in action.option_strings
        }
        self._dest_to_option_strings[action.dest] = action.option_strings
        # if post_type is not None:
        #     self._post_type[action.dest] = post_type
        return action

    def get_arguments_and_config(self):
        args = self.parser.parse_args()
        original_args = copy.copy(args)
        graph = self._default_graph.copy()
        for dest, value in vars(args).items():
            if value is not None:
                graph.nodes[dest]["value"] = value
                graph.nodes[dest]["priority"] = -2
        for v in graph.nodes:
            try:
                if v is not self.default_node and graph.nodes[v]["value"] \
                   is not None:
                    graph.remove_edges_from(list(graph.in_edges(v)))
            except KeyError:
                from IPython import embed; embed()
        propagate_defaults(graph)
        config_path = graph.nodes["input_config"]["value"]
        config = self._config_class()
        if config_path is not None:
            config = self._config_class.yaml_load(config_path)
            for field in self._config_class.__dataclass_fields__:
                value = getattr(config, field)
                if value is not None \
                   and graph.nodes[field]["value"] is None:
                    graph.nodes[field]["value"] = value
                    graph.nodes[field]["priority"] = -1
            propagate_defaults(graph)
        for node in graph.nodes:
            if node is not self.default_node:
                if graph.nodes[node]["required"] \
                   and graph.nodes[node]["value"] is None:
                    err = f"Missing required value {node}."
                    try:
                        err += " You can provide it with {}.".format(
                            self._dest_to_option_strings[node][0]
                        )
                    except KeyError:
                        pass
                    self.parser.error(err)
                if hasattr(args, node):
                    setattr(args, node, graph.nodes[node]["value"])
                if hasattr(config, node):
                    setattr(config, node, graph.nodes[node]["value"])
        sc = None
        try:
            sc = set(args._config_show_config)
        except (AttributeError, TypeError):
            pass
        if sc is not None:
            try:
                show_config_format = args._config_show_config_format
            except AttributeError:
                if sc == {"config"}:
                    show_config_format = "yaml"
                else:
                    show_config_format = "dict"
            to_print = []
            for name in ("original_args", "args"):
                if name in sc:
                    rep = f"Could not represent in format {show_config_format}."
                    if show_config_format == "dict":
                        rep = pprint.pformat(vars(locals()[name]))
                    elif show_config_format == "json":
                        rep = json.dumps(vars(locals()[name]), indent=4)
                    elif show_config_format == "yaml":
                        rep = dump(vars(locals()[name]), Dumper=Dumper)
                    header = "="*(len(name)+2)
                    to_print.append(
                        f"# {header}\n#  {name} \n# {header}\n\n{rep}"
                    )
            if "config" in sc:
                rep = f"Could not represent in format {show_config_format}."
                if show_config_format == "dict":
                    rep = pprint.pformat(dataclasses.asdict(config))
                elif show_config_format == "json":
                    rep = json.dumps(config.marshal(), indent=4)
                elif show_config_format == "yaml":
                    rep = dump(config.marshal(), Dumper=Dumper)
                header = "="*len(" config ")
                to_print.append(f"# {header}\n#  config \n# {header}\n\n{rep}")
            print("\n\n".join(to_print))
            sys.exit(0)
        return original_args, args, config

    def expose_config_field(
            self,
            f,
            aliases=[],
            default=None,
            help=None,
            action=None,
            default_alias=True,
            nargs=None,
            positional=False,
            #use_post_type=None,
            **kwargs
    ):
        # if use_post_type is None:
        #     use_post_type = self._use_post_type
        if default_alias:
            aliases = ["--{}".format(f.replace("_", "-"))] + aliases
        if positional:
            aliases = [f]
        if default is None:
            default = self._config_class.__dataclass_fields__[f].default
        if help is None:
            try:
                help = self._config_class.__dataclass_fields__[
                    f
                ].metadata["description"]
            except KeyError:
                pass
        type_ = str
        try:
            type_ = self._config_class.__dataclass_fields__[f].type
        except AttributeError:
            pass
        try:
            action = self._type_to_action[type_]
        except KeyError:
            pass
        #print(aliases)
        if action is None:
            try:
                kwargs["unmarshal"] = self._config_class.__dataclass_fields__[
                    f
                ].unmarshal
                action = UnmarshallingStoreAction
            except AttributeError:
                pass
        if nargs is None:
            nargs = detect_nargs(type_, self._multi_nargs_type)
        if not positional:
            kwargs["dest"] = f
        if action not in ["store_true", "store_false"]:
            kwargs["nargs"] = nargs
        return self.add_argument(
            *aliases,
            default=default,
            help=help,
            action=action,
            **kwargs
        )

    def set_required(self, field, value=True):
        self._default_graph.nodes[field]["required"] = value

    def expose_fields(self, fields, **kwargs):
        for f in fields:
            self.expose_config_field(f, **kwargs)
        
class RNACliqueConfigArgumentManager(ConfigArgumentManager):
    default_out_dirs = {
        "top_genes_dir": "od1",
        "tables_dir": "od2",
        "cache_dir": "db_cache",
    }

    default_out_files = {
        "graph": "graph.pkl",
        "matrix": "distance_matrix.h5",
        "output_config": "config.yaml",
    }        

    default_in_names = {
        "input_config": "config.yaml"
    }

    default_aliases = {
        "top_genes": ["-n"],
        "top_matches": ["-N"],
        "transcripts_name": ["--transcripts", "-t"],
        "top_genes_dir": ["--out-dir-1", "-O1"],
        "tables_dir": ["--out-dir-2", "-O2"],
        "transcript_id_regex": ["--pattern", "-p"],
        "evalue": ["-e"],
        "graph": ["-g", "--output-graph"],
        "jobs": ["-j"],
        "cache_dir": ["-C"],
        "output_dir": ["-O"],
        "input_dirs": ["-I"],
        "matrix": ["-m"],
        "title": ["-T"],
    }
    
    def __init__(self, *args, default_deps={}, **kwargs):
        super().__init__(
            RNACliqueConfig,
            *args,
            default_deps=default_deps,
            **kwargs
        )
        for name, filename in type(self).default_out_names.items():
            self._add_out_file_deps(name, filename)
        for name, filename in type(self).default_in_names.items():
            self._add_in_file_deps(name, filename)
        self.add_defaults("version", get_version(), -3)
        self.add_defaults(
            "title",
            {
                ("output_dir",): lambda output_dir: output_dir.name
            }
        )

    def _add_in_file_deps(self, name, filename):
        def inner(output_dir):
            f = output_dir / filename
            if f.exists():
                return f
            return None
        self.add_defaults(name, {("output_dir",): inner})

    def _add_out_file_deps(self, name, filename):
        def inner(output_dir):
            return output_dir / filename
        self.add_defaults(name, {("output_dir",): inner})

    def expose_fields_with_default_aliases(self, *fields, **kwargs):
        for f in fields:
            self.expose_config_field(
                f,
                aliases=self.default_aliases.get(f, []),
                **kwargs
            )

    def add_output_config_argument(self, aliases=["--output-config", "-c2"]):
        self.add_argument(
            *aliases,
            dest="output_config",
            type=Path,
        )

    @classmethod
    def make_output_dirs(cls, config):
        try:
            config.output_dir.mkdir(exist_ok=True)
        except AttributeError:
            pass
        for d in cls.default_out_dirs:
            getattr(config, d).mkdir(exist_ok=True)
        

RNACliqueConfigArgumentManager.default_out_names = ChainMap(
    RNACliqueConfigArgumentManager.default_out_dirs,
    RNACliqueConfigArgumentManager.default_out_files
)
