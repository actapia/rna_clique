import functools

import Bio.Phylo
import numpy as np
import pandas as pd
import matplotlib as mpl
import dendropy

from io import StringIO
from numbers import Number
from typing import Callable, Union, Optional, Any, Iterator, TypeVar
from collections.abc import Sequence, Mapping

from matplotlib import pyplot as plt

from .plots import _transform_ax, BasicCompositeTransform
from ..identity import id_

def tril_jagged(mat: np.ndarray) -> list[list[Number]]:
    """Get lower triangle of a matrix as a jagged 2-dimensional list."""
    return [r[:(i+1)].tolist() for (i, r) in enumerate(np.tril(mat))]

T = TypeVar("T")

def draw_tree(
        tree: Bio.Phylo.BaseTree.Tree,
        blank_nonterminals: bool = True,
        clades: Optional[dict[T, Bio.Phylo.BaseTree.Clade]] = None,
        colors: Optional[
            Union[
                mpl.colors.ListedColormap,
                Sequence[tuple[int, int, int]],
                Mapping[T, tuple[int, int, int]]
            ]
        ] = None,
        ax = None
) -> mpl.axes.Axes:
    """Draw a BioPython tree using matplotlib.

    This function builds on the Bio.Phylo.draw function included in BioPython to
    enable coloring clades. The function also allows automatic blanking of
    non-terminals.

    The clades parameter, if provided, should be a dict mapping values (of any
    type) to clades belonging to the given tree. These clades will be colored by
    using the provided colors parameter.

    The colors parameter can provide colors by index or by value. To color by
    index, an mpl.colors.ListedColormap object or a Sequence (e.g., a list) of
    RGB255 (0--255 range) tuples should be provided. In that case, colors will
    be assigned to clades in the same order that they appear in the clades
    dict. To color by value, a Mapping (e.g., a dict) object mapping the values
    associated with the clades to RGB255 colors should be provided.

    Parameters:
        tree:                      The phylogenetic tree to draw.
        blank_nonterminals (bool): Give non-terminals blank labels.
        clades (dict):             dict mapping values to clades
        colors:                    Object for obtaining colors for clade values.
        ax:                        Matplotlib Axes on which to draw.

    Returns:
        The matplotlib Axes object on which the tree was drawn.
    """
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

def bbox_max_x(bbox: mpl.transforms.BboxBase) -> float:
    """Get the maximum x coordinate of a bounding box."""
    return max(bbox.bounds[0] + bbox.bounds[2], bbox.bounds[0])

def bbox_min_y(bbox: mpl.transforms.BboxBase) -> float:
    """Get the minimum y coordinate of a bounding box."""    
    return min(bbox.bounds[1], bbox.bounds[1] + bbox.bounds[3])

def bbox_max_y(bbox: mpl.transforms.BboxBase) -> float:
    """Get the maximum y coordinate of a bounding box."""    
    return max(bbox.bounds[1], bbox.bounds[1] + bbox.bounds[3])

def draw_clade_labels(
        ax: mpl.axes.Axes,
        clades: dict[T, Bio.Phylo.BaseTree.Clade],
        colors: Optional[
            Union[
                mpl.colors.ListedColormap,
                Sequence[tuple[int, int, int]],
                Mapping[str, tuple[int, int, int]]
            ]
        ] = None,
        line_padding: float = 0.036,
        cap_width: float = 0.02,
        text_padding: float = 0.023,
        make_label: Callable[[T], str] = id_
):
    """Draw labels for specified clades at the right of the plot.

    This function draws vertical capped line segments on the right side of a
    plot created with draw_tree. Each such line segment denotes a minimal
    interval on the vertical axis containing all terminals of a given clade.

    The clades to be labeled are specified in a dict provided as a parameter to
    this function. The dict should map values identifying the clades to the
    clades themselves. By default, the label drawn with the line segment for a
    clade is simply the value with which the clade is associated in the clades
    dict. To use a different string as a label, one can provide the make_label
    parameter, which should be a unary function mapping keys from the clade dict
    to the labels (strings) to use.

    The positioning and dimensions of the line segments drawn can be controlled
    with the line_padding, cap_width, and text_padding parameters. All three are
    expected to be provided in axis coordinates.

    The line_padding parameter dictates how far the drawn line segments should
    be offset horizontally from the right edges of the existing plot's terminal
    labels.

    cap_width controls the size of the caps of the line segmeents drawn.

    text_padding controls the horizontal distance from the drawn lines to the
    clade text labels.

    As in draw_tree, colors may be provided to color the clade labels.  The
    colors parameter can provide colors by index or by value. To color by index,
    an mpl.colors.ListedColormap object or a Sequence (e.g., a list) of RGB255
    (0--255 range) tuples should be provided. In that case, colors will be
    assigned to clades in the same order that they appear in the clades dict. To
    color by value, a Mapping (e.g., a dict) object mapping the values
    associated with the clades to RGB255 colors should be provided.

    Parameters:
        ax:
        clades (dict):
        colors:               
        line_padding (float): Left padding for lines relative to existing text.
        cap_width (float):    Size of caps to draw on line segments.
        text_padding (float): Left padding for text labels relative to lines.
        make_label:           Function to make labels from clade dict keys.
    """
    def ax_to_data_x(x: float) -> float:
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
        min_y = min(
            bbox_min_y(bbox_dict[t.name]) for t in clade.get_terminals()
        )
        max_y = max(
            bbox_max_y(bbox_dict[t.name] )for t in clade.get_terminals()
        )
        mid_y = (min_y + max_y)/2
        #print(label, max_x, (min_y, mid_y, max_y))
        ax.plot([max_x, max_x], [min_y, max_y], clip_on=False, color=line_color)
        # Plot caps.
        cap_x = [max_x - cap_width/2, max_x + cap_width/2]
        ax.plot(cap_x, [min_y, min_y], clip_on=False, color=line_color)
        ax.plot(cap_x, [max_y, max_y], clip_on=False, color=line_color)
        ax.text(
            max_x + text_padding, mid_y,
            make_label(label),
            horizontalalignment="left",
            verticalalignment="center",
            color=text_color
        )
    ax.get_figure().tight_layout()

def phylo_to_dendropy(
        tree: Bio.Phylo.BaseTree.Tree,
        ns: Optional[dendropy.TaxonNamespace] = None
) -> dendropy.Tree:
    """Convert a BioPython Tree into a Dendropy Tree.

    Optionally, a TaxonNamespace can be provided to ensure that taxa with
    the same name point to the same Taxon when making multiple trees this
    way. If no such namespace is provided, a new one will be created with the
    new tree.

    Parameters:
        tree: BioPython Tree object to convert.
        ns:   Dendropy TaxonNamespace to use for taxa in the new tree.

    Returns:
        An equivalent Dendropy Tree using any provided namespace.
    """
    sio = StringIO()
    Bio.Phylo.write(tree, sio, format="newick")
    return dendropy.Tree.get(
        data=sio.getvalue(),
        schema="newick",
        taxon_namespace=ns
    )

def multi_phylo_to_dendropy(
        *trees: Bio.Phylo.BaseTree.Tree,
        ns: Optional[dendropy.TaxonNamespace] = None
) -> Iterator[dendropy.Tree]:
    """Convert multiple BioPython Trees to Dendropy Trees using one namespace.

    If no namespace is provided as a parameter to the function, one will be
    created with the first tree converted.

    Parameters:
        ns: The DendroPy namespace to use for the trees' taxa.
    """
    ns = None
    for tree in trees:
        new = phylo_to_dendropy(tree, ns=ns)
        ns = new.taxon_namespace
        yield new

def get_clades(
        tree: Bio.Phylo.BaseTree.Tree,
        sample_metadata: pd.DataFrame,
        sample_name_column: str,
        group_by: str | list[str]
) -> Iterator[tuple[Any, Bio.Phylo.BaseTree.Clade]]:
    """Get maximal clades containing terminals with exactly one metadata value.

    This function tries to find the largest clades such that every terminal of
    each clade has the same value for the given column(s) in the metadata, and
    every terminal with that value belongs to the clade.

    Parameters:
        tree:                        Phylogenetic tree for which to get clades.
        sample_metadata:             DataFrame of clade terminal metadata.
        sample_name_column (str):    Column containing sample names.
        group_by:                    Column(s) from which to get values.
    """
    terminals = {c.name: c for c in tree.get_terminals()}
    for group, df in sample_metadata.groupby(group_by):
        clade = tree.common_ancestor(
            *(terminals[n] for n in df[sample_name_column])
        )
        if clade.count_terminals() == df.shape[0]:
            yield group, clade
