import argparse
import importlib.util
import sys
import functools
import re
import numbers

import pandas as pd

from pathlib import Path
from typing import Optional

from IPython import embed

from .. import config as config_module
from ..app import eprint
from .markdown import MarkdownDocument
from .docs import unoptional, get_type_name, column_to_text

checks = ["missing-argument-description", "missing-program-description"]

def build_parser():
    parser = config_module.ArgumentManager()
    parser.add_argument(
        "modules",
        nargs="+",
        type=Path,
        help="Modules for which to automatically generate usage info."
    )
    parser.add_argument(
        "--depth",
        "-d",
        type=int,
        default=1,
        help="Initial Markdown heading depth.",
    )
    parser.add_argument(
        "--width",
        "-w",
        type=int,
        default=80,
        help="Text fill column (maximum line width).",
    )
    parser.add_argument(
        "--ignore-missing-parsers",
        "-i",
        action="store_true",
        help="Ignore modules that lack a build_parser function.",
    )
    parser.add_argument(
        "--check",
        "-c",
        nargs="*",
        action=config_module.UnmarshallingStoreAction,
        unmarshal=set,
        default=set(),
        const=checks,
        choices=checks,
        help="Checkers to run."
    )
    return parser

def import_file(module):
    # https://stackoverflow.com/a/67692
    spec = importlib.util.spec_from_file_location(module.stem, module)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[module.stem] = mod
    spec.loader.exec_module(mod)
    return mod    

nargs_format = {
    "+": r"$\ge 1$",
    "*": r"$\ge 0$",
    "?": r"$0--1$"
}

sci_float_re = re.compile("(.*)e(.*)")

def format_real(const: numbers.Real, marshal=str) -> str:
    return sci_float_re.sub(r"\1 \\times 10^{\2}", str(marshal(const)))

def format_constant(const, marshal) -> str:
    if isinstance(const, numbers.Number) and not isinstance(const, bool):
        parts = []
        if const.real:
            parts.append(format_real(const.real, marshal))
        if const.imag:
            parts.append(format_real(const.imag, marshal) + "i")
        return "${}$".format(" + ".join(parts))
    return MarkdownDocument.md_code(str(marshal(const)))

def format_description(rule, marshal=str):
    match rule.description.type_:
        case config_module.DescriptionType.plain: 
            return str(rule.description)
        case config_module.DescriptionType.constant:
            try:
                return format_constant(rule.value, marshal)
            except AttributeError:
                try:
                    return format_constant(rule(), marshal)
                except TypeError:
                    pass
    return MarkdownDocument.md_code(str(rule.description))

def get_preferred_option_string(option_strings) -> Optional[str]:
    try:
        return next(s for s in option_strings if s.startswith("--"))
    except StopIteration:
        return option_strings[0]

def summarize_cli_args(parser):
    positional_rows = []
    option_rows = []
    for dest, action in parser._dest_to_action.items():
        row = {}
        if dest in parser._config_class.__dataclass_fields__:
            row["config_option"] = MarkdownDocument.md_code(dest)
            row["sort_name"] = dest
        if action.option_strings:
            row.setdefault(
                "sort_name",
                get_preferred_option_string(action.option_strings)
            )
            for s in action.option_strings:
                if s[1:].startswith("-"):
                    row.setdefault("long_name", MarkdownDocument.md_code(s))
                else:
                    row.setdefault("short_name", MarkdownDocument.md_code(s))
                if "short_name" in row and "long_name" in row:
                    break
            option_rows.append(row)
        else:
            # row["name"] = dest
            positional_rows.append(row)
        if action.nargs is not None:
            row["nargs"] = nargs_format.get(action.nargs, f"${action.nargs}$")
        elif action.const:
            row["nargs"] = "$0--1$"
        else:
            row["nargs"] = "$1$"
        row["nargs"] = str(row["nargs"])
        edge_groups = [
            (k, list(v)) for (k, v) in
            config_module.get_in_edge_groups(
                parser._default_graph,
                dest
            )
        ]
        try:
            marshal = parser._config_class.__dataclass_fields__[
                dest
            ].marshal
        except KeyError:
            marshal = str
        if action.default:
            if action.default != argparse.SUPPRESS:
                row["default"] = format_constant(
                    action.default,
                    marshal=marshal
                )
        elif edge_groups:
            edge_data = parser._default_graph.edges[edge_groups[0][1][0]]
            
            fd = functools.partial(
                format_description,
                edge_data["function"],
                marshal=marshal,
            )
            if len(edge_groups) > 1:
                if edge_data["args"] == () \
                   and edge_data["function"].description is not None:
                    row["default"] = fd()
            else:
                args = []
                for arg in edge_data["args"]:
                    if arg.startswith("_"):
                        try:
                            args.append(
                                get_preferred_option_string(
                                    parser._dest_to_action[arg].option_strings
                                )
                            )
                        except IndexError:
                            pass
                    else:
                        args.append(arg)
                if edge_data["function"].description is not None:
                    row["default"] = fd()
                elif args != ():
                    args_str = text_list(
                        map(
                            MarkdownDocument.md_code,
                            args
                        ),
                        "and"
                    )
                    row["default"] = f"Depends on {args_str}"
            row.setdefault("default", "Dynamic")
        try:
            row["required"] = parser._default_graph.nodes[dest]["required"]
        except KeyError:
            pass
        row["description"] = action.help
        try:
            row["type"] = unoptional(
                parser._config_class.__dataclass_fields__[dest].type
            )
        except KeyError:
            row["type"] = action.type
            if not row["type"] and action.nargs != 0:
                row["type"] = str
            if action.nargs not in (None, "?", 0, 1):
                row["type"] = list[row["type"]]
        if row["type"]:
            type_name = get_type_name(row["type"])
            row["type"] = MarkdownDocument.md_code(type_name)
        row["choices"] = action.choices
        if row["choices"]:
            try:
                row["choices"] = text_list(
                    map(
                        MarkdownDocument.md_code,
                        row["choices"]
                    ),
                    "or"
                )
            except TypeError:
                embed()
        if action.const is not None:
            row["const"] = format_constant(action.const, marshal)
    return (
        pd.DataFrame(positional_rows).rename_axis("position"),
        pd.DataFrame(option_rows),
    )

def text_list(l, end: str) -> str:
    l = list(l)
    if len(l) == 1:
        return str(l[0])
    if len(l) == 2:
        return f" {end} ".join(map(str, l))
    return "{}, {} {}".format(
        ", ".join(map(str, l[:-1])),
        end,
        str(l[-1])
    )

def remove_empty_columns(df):
    return df[df.columns[~df.isnull().all()]]

new_column_names = {
    "nargs": "Argument count",
    "default": "Default value",
    "const": "Default value (flag only)"
}

def default_columns(df, l):
    col = df[l[0]].copy()
    for c in l:
        ix = col.isnull()
        col.loc[ix] = df[c].loc[ix]
    return col

def sort_existing(df, by, *args, **kwargs):
    return df.sort_values([b for b in by if b in df.index], *args, **kwargs)    
    
def main():
    _, args = build_parser().get_arguments()
    md = MarkdownDocument(depth=args.depth)
    with md.section("Command-line usage guide"):
        for module_file in args.modules:
            module = import_file(module_file)
            try:
                bp = module.build_parser
            except AttributeError as e:
                if args.ignore_missing_parsers:
                    eprint(f"Missing parser in {module_file}.")
                    continue
                raise e
            parser = bp()
            with md.section(MarkdownDocument.escape(str(module_file))):
                if parser.parser.description:
                    md.paragraph(parser.parser.description)
                elif "missing-program-description" in args.check:
                    eprint(f"Missing program description for {module_file}.")
                positional_args, optional_args = summarize_cli_args(parser)
                if "missing-argument-description" in args.check:
                    if not positional_args.empty:
                        for ix, row in positional_args.loc[
                                positional_args["description"].isna()
                        ].iterrows():
                            eprint(
                                "Missing description for positional argument "
                                f"{ix} in {module_file}"
                            )
                    if not optional_args.empty:
                        for ix, row in optional_args.loc[
                                optional_args["description"].isna()
                        ].iterrows():
                            eprint(
                                "Missing description for {} in {}".format(
                                    row["sort_name"],
                                    module_file,
                                )
                            )
                     
                #embed()
                optional_args["required"] = optional_args["required"].replace(
                    {
                        False: "No",
                        True: "Yes"
                    }
                )
                if not positional_args.empty:
                    with md.section("Positional arguments"):
                        md.paragraph(
                            positional_args.fillna(
                                ""
                            ).rename_axis(
                                "Position"
                            ).filter(
                                [
                                    "config_option",
                                    "description",                       
                                    "nargs",
                                    "type",
                                ]
                            ).rename(
                                columns=new_column_names
                            ).rename(
                                columns=column_to_text
                            ).to_markdown(),
                            wrap=False
                        )
                if not optional_args.empty:
                    with md.section("Options"):
                        #embed()
                        md.paragraph(
                            sort_existing(
                                optional_args,
                                [
                                    "config_option",
                                    "sort_name"
                                ]
                            ).fillna("").filter(
                                [
                                    "config_option",
                                    "long_name",
                                    "short_name",
                                    "description",
                                    "nargs",
                                    "type",
                                    "choices",
                                    "default",
                                    "const",
                                    "required",
                                ]
                            ).rename(
                                columns=new_column_names
                            ).rename(
                                columns=column_to_text
                            ).to_markdown(index=None),
                            wrap=False
                        )

if __name__ == "__main__":
    main()
