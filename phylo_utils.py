import Bio.Phylo
import numpy as np
import functools
import matplotlib as mpl

from plots import _transform_ax, BasicCompositeTransform

from nice_colorsys import *
from nice_colorsys.rgb255 import rgb255

from numbers import Number
from matplotlib import pyplot as plt
from Bio.Phylo import BaseTree

from typing import Callable, Union, Optional
from collections.abc import Iterable, Sequence, Mapping

def tril_jagged(mat: np.ndarray) -> list[list[Number]]:
    """Get lower triangle of a matrix as a jagged 2-dimensional list."""
    return [r[:(i+1)].tolist() for (i, r) in enumerate(np.tril(mat))]

def draw_tree(
        tree: BaseTree,
        blank_nonterminals: bool = True,
        clades: Optional[dict] = None,
        colors: Optional[
            Union[
                mpl.colors.ListedColormap,
                Sequence[tuple[int, int, int]],
                Mapping[str, tuple[int, int, int]]
            ]
        ] = None,
        ax = None
):
    if blank_nonterminals:
        for c in tree.get_nonterminals():
            c.name = None
    if ax is None:
        ax = plt.axes()
    color_dict = {}
    if colors is not None:
        for i, (label, clade) in enumerate(clades.items()):
            color = None
            try:
                color = colors.colors[i]
            except AttributeError:
                try:
                    color = colors[label]
                except IndexError:
                    color = colors[i]
                except TypeError:
                    pass
            if color is not None:
                #print(color)
                #color = rgb255(*color).as_hex()
                for term in clade.get_terminals():
                    color_dict[term.name] = color
                clade.color = color
    Bio.Phylo.draw(tree, do_show=False, axes=ax, label_colors=color_dict)
    ax.axis("off")
    return ax

def bbox_max_x(bbox):
    return max(bbox.bounds[0] + bbox.bounds[2], bbox.bounds[0])

def bbox_min_y(bbox):
    return min(bbox.bounds[1], bbox.bounds[1] + bbox.bounds[3])

def bbox_max_y(bbox):
    return max(bbox.bounds[1], bbox.bounds[1] + bbox.bounds[3])

def draw_clade_labels(
        ax,
        clades: dict,
        colors: Optional[
            Union[
                mpl.colors.ListedColormap,
                Sequence[tuple[int, int, int]],
                Mapping[str, tuple[int, int, int]]
            ]
        ] = None,
        line_padding: float = 0.036,
        cap_width: float = 0.02,
        text_padding: float = 0.023
):
    def ax_to_data_x(x):
        return ax_to_data([x],0)[0] - ax_to_data([0],0)[0]
    ax_to_data = functools.partial(
        _transform_ax,
        BasicCompositeTransform(ax.transAxes, ax.transData.inverted())
    )
    line_padding = ax_to_data_x(line_padding)
    text_padding = ax_to_data_x(text_padding)
    cap_width = ax_to_data_x(cap_width)
    bbox_dict = {
        t.get_text().strip(): t.get_window_extent(
            ax.get_figure().canvas.get_renderer()
        ).transformed(ax.transData.inverted())
        for t in ax.texts
    }
    max_x = max(bbox_max_x(b) for b in bbox_dict.values()) + line_padding
    for i, (label, clade) in enumerate(clades.items()):
        line_color = "0.8"
        text_color = "black"
        color = None
        try:
            color = colors.colors[i]
        except AttributeError:
            try:
                color = colors[label]
            except IndexError:
                color = colors[i]
            except TypeError:
                pass
        if color is not None:
            line_color = color
            text_color = color
        #max_x = max(pos_dict[t.name][0] for t in clade.get_terminals())
        min_y = min(bbox_min_y(bbox_dict[t.name]) for t in clade.get_terminals())
        max_y = max(bbox_max_y(bbox_dict[t.name] )for t in clade.get_terminals())
        mid_y = (min_y + max_y)/2
        #print(label, max_x, (min_y, mid_y, max_y))
        ax.plot([max_x, max_x], [min_y, max_y], clip_on=False, color=line_color)
        # Plot caps.
        cap_x = [max_x - cap_width/2, max_x + cap_width/2]
        ax.plot(cap_x, [min_y, min_y], clip_on=False, color=line_color)
        ax.plot(cap_x, [max_y, max_y], clip_on=False, color=line_color)
        ax.text(
            max_x + text_padding, mid_y,
            label,
            horizontalalignment="left",
            verticalalignment="center",
            color=text_color
        )
