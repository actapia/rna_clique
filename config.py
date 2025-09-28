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
from typing import Optional, Any
from collections import deque, ChainMap
from collections.abc import Collection, Iterable
from pathlib import Path

from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from marshalling_dataclass import (
    marshalling_dataclass,
    marshalling_field,
    MarshallingDataclassBase,
)
from transcripts import default_gene_re

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
    """Dataclass containing configuration data for RNA-clique programs."""
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
        """Load configuration from a YAML file.

        Parameters:
            path: Path to the YAML file from which to load the configuration.
        """
        with open(path, "r") as f:
            return cls.from_marshalled_representation(
                load(
                    f,
                    Loader=Loader
                )
            )

    def mark_finish(self):
        """Set the finished attribute to the current date and time."""
        self.finished = datetime.now()

    def yaml_save(self, path: Path):
        """Save the configuration in a YAML file.

        Parameters:
            path: The path to the YAML file to which to save the configuration.
        """
        with open(path, "w") as f:
            f.write(dump(self.marshal(hide_none=True), Dumper=Dumper))

    @classmethod
    def from_args(cls, args):
        """Create a configuration from an argparse Namespace.

        If there is an input_config attribute for the namespace, and it is not
        None, the configuration file will be loaded from the specified path. Any
        remaining attributes of the namespace that are not None and that match
        fields of the configuration dataclass will be used to override values in
        the loaded configuration.

        Parameters:
            args: The argparse namespace on which to base the configuration.
        """
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

def detect_nargs(t, multi_type: str="+") -> Optional[str]:
    """Detect whether an attribute type should correspond to multiple arguments.

    This function determines whether an attribute with a given type annotation
    should correspond to a command-line option that accepts multiple arguments
    via argparse's nargs keyword argument.

    Parameters:
        t:                The type annotation for the attribute.
        multi_type (str): nargs value to use when accepting multiple arguments.

    Returns:
        The value to use for argparse nargs for the attribute type.
    """
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

def get_edge_priority(graph: nx.Graph, edge) -> int:
    """Get the priority of an edge in the given graph."""
    return graph.edges[edge]["priority"]

def update_node(graph: nx.Graph, origin) -> bool:
    """Update a configuration graph node with higher priority values.

    This function tries to update the value of a node in the configuration
    graph. It checks whether there is a group of edges incident to the node that
    can be used to update the node to a value with priority lower than the
    current priority.

    Parameters:
        graph:  The configuration graph in which to perform the update.
        origin: The node to update.

    Returns:
        A boolean indicating whether the node was updated.
    """
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
    

def propagate_defaults(graph: nx.Graph):
    """Update node values in the graph to minimize priority.

    This function propagates values through the edges of the graph in such a way
    that each node is assigned a value with minimal possible priority.

    Parameters:
        graph: The configuration graph in which to update all node values.
    """
    con = nx.condensation(graph)    
    for scc in nx.topological_sort(con):
        st = deque()
        for pred in con.predecessors(scc):
            # print("pred is", con.nodes[pred]["members"])
            # print("scc is", con.nodes[scc]["members"])
            for node in nx.node_boundary(
                    graph,
                    con.nodes[pred]["members"],
                    con.nodes[scc]["members"]
            ):
                
                if update_node(graph, node):
                    st.append(node)
        subgraph = graph.subgraph(con.nodes[scc]["members"])                
        # st = deque(
        #     (n, True)
        #     for n in subgraph.nodes if subgraph.nodes[n]["value"] is not None
        # )
        while st:
            origin = st.pop() #, changed = st.pop()
            changed = update_node(subgraph, origin)
            if changed:
                for edge in subgraph.out_edges(origin):
                    dest = edge[1]
                    if subgraph.edges[edge]["priority"] <= \
                       subgraph.nodes[dest]["priority"]:
                        st.append(dest)

class DefaultDiGraph(nx.DiGraph):
    """Directed graph for configuration (defaults) graph.

    Unlike nx.DiGraph, this subclass accepts default attributes to assign to
    nodes and edges.
    """
    def __init__(self, *args, node_defaults={}, edge_defaults={}, **kwargs):
        super().__init__(*args, **kwargs)
        self.node_defaults = node_defaults
        self.edge_defaults = edge_defaults

    def node_attr_dict_factory(self) -> dict:
        return dict(self.node_defaults)

    def edge_attr_dict_factory(self) -> dict:
        return dict(self.edge_defaults)

class MultiArgConstAction(argparse.Action):
    """Argparse action that allows use of const with nargs='*' or nargs='+'"""
    def _get_values(self, values):
        if not values:
            return self.const
        else:
            return values
        
    def __call__(self, parser, namespace, values, option_string):
        setattr(namespace, self.dest, self._get_values(values))

class UnmarshallingStoreAction(MultiArgConstAction):
    """Argparse action that performs unmarshalling."""
    # TODO: Remind myself why this is necessary instead of using the type kwarg.
    def __init__(self, *args, unmarshal=lambda x: x, **kwargs):        
        super().__init__(*args, **kwargs)
        self.unmarshal = unmarshal

    def _get_values(self, values):
        return self.unmarshal(super()._get_values(values))


class ConfigArgumentManager[T: MarshallingDataclassBase]:
    """Class that manages configuration via config files and CLI arguments.

    This class manages configuration controlled by configuration files and
    command-line arguments and handles the interplay between the two. The
    class's interface is similar to argparse's ArgumentParser, but this class
    can also load settings from configuration files and can automatically expose
    settings as command-line arguments in a sensible way using the configuration
    dataclass's field information.

    This class also provides a powerful system for inferring setting values from
    other settings. Any configuration variable can be set automatically using
    functions that accept other variables' values as arguments. Additionally,
    multiple such functions can be provided with different priorities so that
    certain ways of automatically assigning the value are preferred.

    Arguments accepted by a program need not all correspond to attributes of the
    dataclass used. Arguments not belonging to the config will be returned by
    the get_arguments_and_config method as well.

    Internally, ConfigArgumentManager keeps a directed graph (the "configuration
    graph" or "defaults graph") where vertices represent setting (variable)
    names. An set of edges F from vertices X to a vertex b exist when the value
    for b can be derived automatically when all values of X are known. Each edge
    in F is labelled with the same function that accepts the values of the
    vertices in X as arguments and returns the derived value of b. Each edge (or
    set of edges) has an associated priority; edges with lower priority values
    are preferred when deriving a vertex's value.
    """
    
    # default_node is a "dummy" value that allows configuration graph nodes to
    # inherit constant default values.
    default_node = object()

    # Default actions to use for variables with certain type annotations.
    _type_to_action = {
        bool: "store_true",
        Optional[bool]: "store_true"
    }
    
    def __init__(
            self,
            config_class: type[T],
            *args,
            default_deps: dict[str, Any] = {},
            accept_input: bool = True,
            show_config: bool = True,
            show_config_format: bool = True,
            type_to_action: Optional[dict[Any, str]] = None,
            multi_nargs_type: str = "+",
            **kwargs
    ):
        """Construct a ConfigArgumentManager with the given settings.

        The provided config_class must be a marshalling dataclass; its fields
        will be used to automatically create vertices in the configuration
        graph, though not all fields need to be used.

        Initial dependencies to add can be specified in the default_deps
        argument. Each keys must be the name of a configuration setting for
        which to add defaults, and each value should be either the (independent)
        default value for that setting or another dict. If the value is another
        dict, the keys of that dict should be tuples of other settings that can
        be in automatically assigning the first setting, and the values should
        be the functions that accept those settings as keyword arguments.

        ConfigArgumentManager uses UnmarshallingStoreActions for creating
        argparse arguments by default. These behave similarily to the default
        _StoreAction of argparse. Different default actions to use for specific
        type annotations can be specified using the type_to_action argument of
        this constructor. The actual value for this parameter should be a dict
        mapping type annotations to strings representing actions to take. If no
        value is provided for this keyword argument, the value will default to
        a copy of the _type_to_action dict belonging to the class itself, i.e.,
        ConfigArgumentManager._type_to_action.

        ConfigArgumentManager can automatically add some common arguments whose
        usefulness is independent ofthe specific configuration dataclass
        used. These are enabled by default but can be disabled using the boolean
        arguments to this constructor.

        Parameters:
            config_class (type):       Dataclass to use for storing config.
            default_deps (dict):       Dependencies to add to defaults graph.
            accept_input (bool):       Add argument to accept input config file.
            show_config (bool):        Add argument to show processed config.
            show_config_format (bool): Add argument to set show config format.
            type_to_action (dict):     Mapping from type annotations to actions.
            multi_nargs_type (str):    argparse nargs for multi-arg options.

        """
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

    def add_config_argument(
            self,
            aliases: Iterable[str] = ["--input-config", "-c", "-c1"]
    ):
        """Add a command-line argument for providing an input config file.
        
        Parameters:
            aliases: Option strings to use for input config argument.
        """
        self.add_argument(
            *aliases,
            dest="input_config",
            type=Path
        )

    def add_show_config_argument(
            self,
            aliases: Iterable[str] = ["--show-config"]
    ):
        """Add a command-line argument for showing the computed configuration.

        The added argument allows the user to show the original parsed
        command-line arguments (immediately after parsing by argparse), the
        updated arguments (after combining with the config file sesttings),
        and/or the configuration (after combining with the command-line
        arguments). By default, only the configuration is shown.

        Parameters:
            aliases: Option strings to use for show config argument.
        """
        self.add_argument(
            *aliases,
            dest="_config_show_config",
            nargs="*",
            action=MultiArgConstAction,
            const=["config"],
            choices=["original_args", "args", "config"]
        )

    @classmethod
    def select_default_config_format(
            cls,
            _config_show_config: list[str]
    ) -> str:
        """Pick a good serializiation for displaying a set of arguments/configs.

        The config should always be serializable as YAML. Since YAML is likely
        the easiest format to read, it is the default format when only the
        config needs to be shown.

        Since the original_args and args might not be cleanly serializable as
        YAML, the default is to use the Python dict repr when viewing those.

        Parameters:
            _config_show_config (list): Arguments/configs to show.

        Returns:
            A good default format for showing the arguments/configs.
        """
        if _config_show_config == ["config"]:
            return "yaml"
        else:
            return "dict"

    def add_show_config_format_argument(
            self,
            aliases: Iterable[str] = ["--show-config-format"]
    ):
        """Add a command-line argument for selecting config/argument format.

        This method add a command-line argument that allows the user to select
        the serialization format when showing the arguments/config. The
        currently supported formats are dict (Python dict repr), json, and yaml.

        Parameters:
            aliases: Option strings to use for the argument.
        """
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

    # def set_priority(
    #         self,
    #         option: str,
    #         pred: tuple[str,...],
    #         priority: int
    # ):
        

    def add_defaults(
            self,
            option: str,
            default,
            next_priority: Optional[int] = None
    ):
        """Add rules (defaults) for constructing a setting automatically.

        The value of the default argument should either be the actual default
        value for the setting or a dict specifying rules for constructing the
        setting from other settings. If a dict is specified, each item
        represents a rule. It should map a tuple of settings to a function that
        accepts those settings as parameters and returns the derived value.

        When a dict is specified for the default, the order of the items in the
        dict matters. Earlier rules have a lower priority value; they take
        precedence over later rules. By default, the first rule is given
        priority one greater than the highest priority rule added so far. If no
        rules have been added prior to calling this method, the first rule is
        given priority 0. Optionally, the priority for the first rule can be
        specified manually with the next_priority argument. Each Subsequent rule
        in the default dict is given priority one greater than the last.

        Unlike the set_defaults method, this method does not remove existing
        defaults.

        Parameters:
            option (str):        Option string of setting to add defaults for.
            default:             Defaults to add for the setting.
            next_priority (int): Priority to assign to first added default.
        """
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

    def set_defaults(
            self,
            option: str,
            default
    ):
        """Set rules (defaults) for constructing a setting automatically.

        The value of the default argument should either be the actual default
        value for the setting or a dict specifying rules for constructing the
        setting from other settings. If a dict is specified, each item
        represents a rule. It should map a tuple of settings to a function that
        accepts those settings as parameters and returns the derived value.

        When a dict is specified for the default, the order of the items in the
        dict matters. Earlier rules have a lower priority value; they take
        precedence over later rules. By default, the first rule is given
        priority one greater than the highest priority rule added so far. If no
        rules have been added prior to calling this method, the first rule is
        given priority 0.

        Unlike the add_defaults method, this method removes any existing
        defaults.

        Parameters:
            option (str): Option string of setting for which to set defaults.
            default:      Defaults to set for the setting.
        """        
        dest = self._option_string_to_dest.get(option, option)
        self._default_graph.remove_edges_from(
            list(self._default_graph.in_edges(dest))
        )
        self.add_defaults(option, default)

    def add_argument(
            self,
            *args,
            required: bool = False,
            default = None,
            #post_type=None,
            **kwargs
    ) -> argparse.Action:
        """Add a command-line argument to be parsed.

        This function operates similarly to the add_argument method of
        argparse.ArgumentParser, and although most arguments/keyword arguments
        are passed as is to argparse, a few are treated specially.

        The required argument is never passed to argparse. Instead,
        ConfigArgumentManager checks if required settings have been specified
        after parsing and propagating defaults through the configuration graph.

        Likewise, the default values (explicit or computed via rules) are
        handled by ConfigArgumentManager and are never passed to argparse.

        A default value can be specified using the default keyword argument to
        this method. The actual value of the default parameter should be either
        a value to use as the default or a dict containing rules for deriving
        the value of the setting from other values. (The former case matches the
        behavior of the default keyword argument in argparse.) If a dict is
        provided, each item of the dict represents a rule. Each key should be a
        tuple of strings naming the settings from which the default value for
        the added setting can be derived, and the corresponding value should be
        a function that accepts those settings as named keyword arguments and
        returns the derived value.

        Parameters:
            required (bool): Whether the setting must always be set.
            default:         Default value or rules for deriving setting value.

        Returns:
            The argparse Action added to the parser.
        """
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

    def get_arguments_and_config(self) -> tuple[
            argparse.Namespace,
            argparse.Namespace,
            T
    ]:
        """Parse and process arguments and config file, if any.

        This function returns three values. The first is the "raw" argparse
        Namespace returned by the parse_args method of the
        argparse.ArgumentParser object maintained by this
        ConfigArgumentManager. The second is the processed argparse Namespace
        that possibly incorporates values for settings derived from other
        settings provided via the command line or a config file. The third and
        final value is the configuration object, which can likewise contain both
        values provided directly in the config file or values derived from the
        command-line arguments and config file settings.

        If any settings marked as required are None in the final configuration
        graph (or, equivalently, are None in both the processed argparse
        Namespace and configuration object), then this method raises an error
        via the ConfigArgumentManager's argparse.ArgumentParser object. If there
        is a command-line argument associated with the first missing setting,
        the first option string for that argument is provided to the user in the
        error message.

        If the user specified the --help option at the command line, the process
        will display argparse's generated help message and terminate before this
        method can return. Likewise, if the --show-config argument is enabled
        and has been specified by the command-line user, this method will print
        the contents of the configuration object (or argparse Namespaces, if
        specified) and terminate.

        A slightly simplified explanation of the operation of this method
        follows.  This function operates by first parsing the command-line
        arguments using argparse. Values specified explicitly in the
        command-line arguments are assigned to the "value" attributes of the
        corresponding settings in a new copy of the configuration graph. These
        values are associated with priority -2, ensuring that they take
        precendence over any values specified in configuration files or derived
        from other settings.

        After the command-line arguments are added to the graph, defaults are
        propagated through the configration graph. This first propagation is
        performed so that actions dependendent on derived settings (e.g.,
        loading the config file) can be performed.

        After the first propagation, the input configuration file, if specified,
        is loaded, and the settings in the configuration file are assigned to
        the value attributes of the corresponding vertices in the configuration
        graph. These values are given priority -1, taking precedence over
        derived values but not command-line settings.

        Finally, propagation is performed again to get the final configuration
        graph. Settings in the final graph matching attributes of the argparse
        Namespace or configuration object overwrite those attributes.

        Returns:
            Unprocessed and processed parsed arguments and the configuration.
        """
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
            f: str,
            aliases: list[str] = [],
            default = None,
            help: Optional[str] = None,
            action: Optional[str] = None,
            default_alias: bool = True,
            nargs: str = None,
            positional: bool = False,
            #use_post_type=None,
            **kwargs
    ) -> argparse.Action:
        """Add a command-line argument for the specified config setting.

        The first argument provided must be the name of an attribute of the
        dataclass used for this ConfigArgumentManager. For arguments that don't
        correspond to any configuration attributes, use the add_argument method
        instead.

        By default, an option string is constructed automatically from the
        attribute name and used in addition to any provided aliases. To only use
        the option strings specified in aliases, the default_alias argument
        should be set to False.

        If no help string is specified in the call to this method, the help
        string will default to the value associated with the "description"
        field's metadata.

        If no action is specified, this method will attempt to select an
        appropriate action by retrieving the value in the type_to_action dict
        (specified in the constructor, defaulting to the class's _type_to_action
        dict) that is associated with the field's
        type annotation. If the type annotation for the field is not present in
        the type_to_action dict, an UnmarshallingStoreAction will be used
        instead.

        If nargs is not specified, this method will try to automatically select
        the number of argv arguments to be accepted by the argparse argument
        using the detect_nargs function. If it is determined that the the
        argparse argument should accept multiple argv arguments, then nargs is
        set to the multi_nargs_type specified in the constructor ("+", by
        default). Otherwise, nargs is None (the default for argparse).

        Although all config settings are exposed as named options by default
        (i.e., arguments indiciated by an option string beginning with - or --
        and followed by the value), it is sometimes more natural to expose
        settings as positional arguments. This can be done by passing True for
        position argument to this method. When a setting is exposed as a
        positional argument, any option strings (default or otherwise) are
        ignored.

        Defaults/rules can be added as a configuration setting is exposed. As in
        add_defaults, the value for the default parameter to this method can
        either be the default value (of any type) or a dict containing rules for
        automatically deriving the setting's value from other settings. If the
        value provided is a dict, each item should map a tuple of strings naming
        settings to a function that accepts values for those settings as keyword
        arguments and returns the derived value for the exposed setting.

        Parameters:
            f (str):              Config field to expose.
            aliases (list):       Additional option strings for exposed field.
            default:              Defaults/rules to use when field unspecified.
            help (str):           Help text for exposed command-line argument.
            action (str):         argparse action to use for exposed argument.
            default_alias (bool): Construct and use a default option string.
            nargs (str):          Number of argv arguments to accept.
            positional (bool):    Whether exposed argument should be positional.

        Returns:
            The created argparse Action for the exposed setting.
        """
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

    def set_required(self, field: str, value: bool = True):
        """Set whether a config setting is required for the program.

        Parameters:
            field (str):  Name of setting for which to set necessity.
            value (bool): Whether the field should be required.
        """
        self._default_graph.nodes[field]["required"] = value

    def expose_fields(self, fields, **kwargs):
        """Expose multiple fields at once using the same parameters."""
        for f in fields:
            self.expose_config_field(f, **kwargs)
        
class RNACliqueConfigArgumentManager(ConfigArgumentManager[RNACliqueConfig]):
    """Manager of configuration and arguments for RNA-clique programs.

    Unlike ConfigArgumentManager, this class contains some attributes and
    methods that are mainly useful for RNA-clique only.
    """
    # Default subdirectories of the output directory.
    default_out_dirs = {
        "top_genes_dir": "od1",
        "tables_dir": "od2",
        "cache_dir": "db_cache",
    }

    # Default output files to be located in the output directory.
    default_out_files = {
        "graph": "graph.pkl",
        "matrix": "distance_matrix.h5",
        "output_config": "config.yaml",
    }        

    # Default input files to be located in the output directory.
    default_in_names = {
        "input_config": "config.yaml"
    }

    # Default aliases for arguments exposing certain configruration fields.
    # These are used in addition to the default option strings created from the
    # field names.
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
    
    def __init__(self, *args, **kwargs):
        """Construct an RNACliqueArgumentManager with the given settings.

        This constructor accepts the same arguments as the parent
        ConfigArgumentManager constructor. It differs from the parent
        constructor in that it automatically adds rules in the configuration
        graph for constructing the certain output and input file/directory paths
        from the root output_directory path.

        This constructor also defaults the version setting to the return value
        of get_version(). This value is given priority -3, so it cannot be
        overridden by command-line arguments or configuration files by default.

        Furthermore, this constructor adds a rule for deriving the "title"
        setting as the name of the output directory.
        """
        super().__init__(RNACliqueConfig, *args, **kwargs)
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

    def _add_in_file_deps(self, name: str, filename: Path | str):
        """Add a rule for deriving an input file setting from the output_dir.

        Parameters:
            name (str): Name of the input file setting.
            filename:   Default filename for the input file.
        """
        def inner(output_dir):
            f = output_dir / filename
            if f.exists():
                return f
            return None
        self.add_defaults(name, {("output_dir",): inner})

    def _add_out_file_deps(self, name: str, filename: Path | str):
        """Add a rule for deriving an output file setting from the output_dir.

        Parameters:
            name (str): Name of the output file setting.
            filename:   Default filename for the output file.
        """
        def inner(output_dir):
            return output_dir / filename
        self.add_defaults(name, {("output_dir",): inner})

    def expose_fields_with_default_aliases(self, *fields: str, **kwargs):
        """Expose multiple fields using the default aliases.

        Additional settings can be provided and apply to all exposed fields.
        """
        for f in fields:
            self.expose_config_field(
                f,
                aliases=self.default_aliases.get(f, []),
                **kwargs
            )

    def add_output_config_argument(
            self,
            aliases: Iterable[str] = ["--output-config", "-c2"]
    ):
        """Add a command-line argument for providing an output config file.
        
        Parameters:
            aliases: Option strings to use for output config argument.
        """        
        self.add_argument(
            *aliases,
            dest="output_config",
            type=Path,
        )

    @classmethod
    def make_output_dirs(cls, config: RNACliqueConfig):
        """Create the output directories if they don't exist already.

        Parameters:
            config: The configuration containing the output directory paths.
        """
        try:
            config.output_dir.mkdir(exist_ok=True)
        except AttributeError:
            pass
        for d in cls.default_out_dirs:
            getattr(config, d).mkdir(exist_ok=True)
        

# This ChainMap map all of the output file and directory fields to their
# corresponding default filenames within the output directory.
RNACliqueConfigArgumentManager.default_out_names = ChainMap(
    RNACliqueConfigArgumentManager.default_out_dirs,
    RNACliqueConfigArgumentManager.default_out_files
)
