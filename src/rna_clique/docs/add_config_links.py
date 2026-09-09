import re
import io
import sys
import functools
import contextlib

import mistletoe
import mistletoe.base_renderer
# import mistletoe.markdown_renderer

import mdpd

from pathlib import Path

from .markdown import MarkdownDocument
from .. import config as config_module

# class MarkdownSectionElement(lxml.etree.ElementBase):
#     def __init__(self, section_name, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#         self.name = section_name
#         self.markdown_elements = []

# def get_section_etree(ast):
#     st = [(MarkdownSectionElement(None), 0)]
#     for x in ast:
#         if isinstance(x, mistletoe.block_elements.Heading):
#             if doc.level ==

# class PlaintextRenderer(mistletoe.base_renderer.BaseRenderer):
#     def render_strong(

def split_target_section(t, delim="/"):
    # rem_re = re.compile(r"^(.*)(\\\\)*(?<!\\)/" + delim)
    # while (m := rem_re.match(t)):
    #     yield match
    strio = io.StringIO()
    accept = True
    for c in t:
        if c == "\\":
            accept = not accept
            if accept:
                strio.write("\\")
        elif c == delim and accept:
            yield strio.getvalue()
            strio = io.StringIO()
        else:
            strio.write(c)
            accept = True
    yield strio.getvalue()

def get_content(elem):
    with mistletoe.base_renderer.BaseRenderer() as renderer:
        return renderer.render(elem).rstrip()

def get_table_line_span(table):
    init = table.line_number
    end = init
    if table.header:
        end = table.header.line_number
    if table.children:
        end = table.children[-1].line_number
    return init - 1, end

heading_re = re.compile("^#* ")

# def get_heading_text(e):
#     return heading_re.sub("", get_content(e))

def get_sections(elems, path, level=None):
    #print(len(elems), path, level)
    if path[0] == "" and level is None:
        #print("A")
        yield from get_sections(elems, path[1:], level=1)
        return
    for i, e in enumerate(elems):
        if isinstance(e, mistletoe.block_token.Heading):
            if level is not None and e.level < level:
                return
            content = get_content(e)
            if (level is None or e.level == level) and \
               (path[0] == "" or content == path[0]):
                #print("Level", level)
                #print("Content", content)
                if len(path) == 1:
                    #print("B")
                    yield e
                else:
                    #print("C")
                    yield from get_sections(
                        elems[i+1:],
                        path[1:],
                        e.level + 1
                    )

def build_parser():
    parser = config_module.ArgumentManager()
    parser.add_argument("source", type=Path)
    parser.add_argument("target", type=Path)
    parser.add_argument("--target-section", default="Settings/")
    parser.add_argument("--column", default="Config option")
    parser.add_argument("-i", "--inplace", nargs="?", const=True)
    return parser

def try_link(available, text):
    try:
        return MarkdownDocument.link(
            text,
            available[get_content(mistletoe.Document(text))]
        )
    except KeyError:
        return text

heading_to_anchor_hyphen = re.compile("[ ]")
heading_to_anchor_empty = re.compile("[.]")

def heading_to_anchor(text):
    return heading_to_anchor_hyphen.sub(
        "-",
        heading_to_anchor_empty.sub(
            "",
            text
        )
    )

def main():
    _, args, config = build_parser().get_arguments_and_config()
    target_section = list(split_target_section(args.target_section))
    with open(args.target, "r") as target_file:        
        sections = [
            get_content(e) for e in get_sections(
                mistletoe.Document(target_file).children,
                target_section
            )
        ]
    source_to_target = args.target.relative_to(args.source.parent)
    section_links = {
        s: "{}#{}".format(
            source_to_target,
            heading_to_anchor(s)
        ) for s in sections
    }
    with open(args.source, "r") as source_file:
        source = source_file.read()
    source_lines = source.splitlines()
    source_doc = mistletoe.Document(source)
    tables = [
        x for x in source_doc.children
        if isinstance(x, mistletoe.block_token.Table)
    ]
    column_to_links = functools.partial(try_link, section_links)
    for table in reversed(tables):
        start, end = get_table_line_span(table)
        df = mdpd.from_md("\n".join(source_lines[start:end]))
        try:
            df[args.column] = df[args.column].apply(column_to_links)
        except KeyError:
            pass
        source_lines[start:end] = df.to_markdown(index=False).splitlines()
    with contextlib.ExitStack() as stack:
        if not args.inplace:
            out = sys.stdout
        else:
            if args.inplace is not True:
                new_source = Path(args.source)
                new_source.name = new_source.name + args.inplace
                args.source.rename(new_source)
            out = open(args.source, "w")                
            stack.enter_context(out)
        for line in source_lines:
            print(line, file=out)    
        
if __name__ == "__main__":
    main()
