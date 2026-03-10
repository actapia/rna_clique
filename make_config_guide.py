import typing

import pandas as pd

import config as config_module

from pathlib import Path

from markdown import MarkdownDocument
from docs import unoptional, get_type_name
from identity import id_

def build_parser():
    parser = config_module.ArgumentManager()
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
        "--config-template",
        "-t",
        type=Path,
    )
    return parser

python_struct_to_yaml = {
    list: "Sequence",
    dict: "Mapping"
}

def python_to_yaml_type(t):
    t = unoptional(t)
    fields = typing.get_args(t)
    type_ = typing.get_origin(t)
    yaml_type = python_struct_to_yaml.get(type_, "Scalar")
    yaml_arg = ""
    yaml_args = []
    yaml_fields = []
    for f in fields:
        try:
            yaml_fields.append(python_to_yaml_type(f))
        except TypeError:
            yaml_fields.append(None)            
    if type_ == list:
        yaml_args = ["of {}"]
    elif type_ == dict:
        yaml_args = ["from {}", "to {}"]
    yaml_arg = " ".join(
        a.format(f) for (a, f) in zip(yaml_args, yaml_fields) if f
    )
    return " ".join(f for f in [yaml_type, yaml_arg] if f)
    

def summarize_config_format(config_class):
    rows = []
    for name, field in config_class.__dataclass_fields__.items():
        row = {}
        row["setting"] = name
        row["python_type"] = MarkdownDocument.md_code(
            get_type_name(
                unoptional(
                    field.type
                )
            )
        )
        if field._marshal == id_:
            yaml_type = field.type
        else:
            yaml_type = field._marshal        
        row["yaml_type"] = python_to_yaml_type(yaml_type)
        row["description"] = field.metadata["description"]
        rows.append(row)
    return pd.DataFrame(rows)

def main():
    _, args = build_parser().get_arguments()
    md = MarkdownDocument(depth=args.depth)
    with md.section("Configuration files"):
        with md.section("Settings"):
            md.paragraph(
                summarize_config_format(
                    config_module.RNACliqueConfig
                ).to_markdown(index=False),
                wrap=False
            )
        try:
            with open(args.config_template, "r") as temp:
                with md.section("Template"):
                    md.code_block(temp.read().strip(), "yaml")
        except TypeError:
            pass

if __name__ == "__main__":
    main()
