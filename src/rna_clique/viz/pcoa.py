import json
import functools

import skbio as skb
import pandas as pd
import matplotlib as mpl

from collections import defaultdict
from typing import Callable, Union, Optional, Any
from collections.abc import Iterable, Sequence, Mapping

from matplotlib import pyplot as plt
from adjustText import adjust_text

from .plots import default_group_label_maker, _keyed_multi_sort, as_tuple
from .confidence_ellipsoid import Ellipsoid, draw_ellipsoid, conf_ellipsoid
from ..identity import id_

def empty_dict(*args, **kwargs) -> dict:
    """Return an empty dict, ignoring any arguments."""
    return {}

try:
    from emperor import Emperor
except ImportError:
    pass

default_marker_style = {"color": "black", "ls": "", "marker": "d"}
scatter_to_line2d = {
    "edgecolors": "mec",
    "linewidth": "mew",
    "s": "ms"
}
scatter_ignore = {"depthshade"}

dim_to_proj = {2: "rectilinear", 3: "3d"}
dimension_default_kwargs = {2: {}, 3: {"depthshade": False}}

def draw_pcoa(
        dis_df: pd.DataFrame,
        sample_metadata: pd.DataFrame,
        group_by: str | Iterable[str],
        sample_name_column: str = "name",
        make_group_label: Callable[[Iterable], str] = default_group_label_maker,
        labelers: Optional[list[Callable[[Any], Any]]] = None,
        colors: Optional[
            Union[
                mpl.colors.ListedColormap,
                Sequence[tuple[int, int, int]],
                Mapping[tuple, tuple[int, int, int]]
            ]
        ] = None,
        legend: bool = True,
        index_to_kwargs: Callable[[int], dict] = empty_dict,
        group_to_kwargs: Callable[
            [str], dict
        ] | Callable[[Iterable[str]], dict] = empty_dict,
        index_group_to_kwargs: Callable[
            [int, str], dict
        ] | Callable[[int, Iterable[str]], dict] = empty_dict,
        order_by: str | Iterable[str] = None,
        sort_key: Optional[
            Union[
                Callable[[pd.Series], pd.Series],
                Iterable[Callable[[pd.Series], pd.Series]]
            ]
        ] = None,
        ellipsoids: bool | Iterable[str] = False,
        make_ellipsoid: Callable[
            [pd.DataFrame],
            Ellipsoid
        ] = conf_ellipsoid(0.95),
        ellipsoid_kwargs: dict[str, Any] = None,
        contribution: bool = True,
        annotate: bool = False,
        adjust: bool = True,
        legend_factors: Optional[list[tuple] | bool] = None,
        default_legend_marker: dict = default_marker_style,
        dropna: bool = True,
        dimensions: int = 2,
        ax: Optional[mpl.axes.Axes] = None,
        #coordinate_sort: bool = False,
        axis_label_kwargs: Optional[dict] = None,
        **scatter_kwargs
) -> skb.stats.ordination.OrdinationResults:
    """Draw a PCoA plot and return the PCoA results.

    This function draws Principal Coordinates Analysis (PCoA) plots. PCoA (also
    called multidimensional scaling or MDS) uses a distance matrix to select
    coordinates for each sample in some fixed number of dimensions so that the
    coordinates best preserve the distances in the distance matrix. A PCoA plot
    plots the coordinates from a PCoA analysis as a scatter plot.

    This function uses scikit-bio's skbio.stats.ordination.pcoa function to
    perform PCoA and uses matplotlib's scatter function to make the scatter
    plot. This function also provides some additional analysis and drawing
    capabilities beyond what is provided by scikit-bio and matplotlib alone.

    Unlike draw_heatmap, this function requires sample metadata because sample
    metadata was judged necessary to make a meaningful PCoA. Nevertheless, there
    may be some instances where sample metadata is not needed, so future
    versions of this function might make sample metadata optional.

    The sample metadata should be a Pandas DataFrame in which each row
    represents a sample, and the columns of the row provide some (arbitrary)
    information about the sample.

    draw_pcoa expects that the rows and columns of the provided distance matrix
    are labeled with some string ID indicating to which sample the row or column
    corresponds. This function should be able to look these IDs up in the sample
    metadata table via some column that provides the same data. Which column
    of the sample metadata provides the sample IDs or names can be specified via
    the sample_name_column parameter.

    Samples can be grouped by sample metadata columns; the set of all samples
    with the same values for those columns forms a group. Although such groups
    have no effect on the coordinates of the samples, they can be used to encode
    information from the sample metadata visually. The columns of the sample
    metadata to use for grouping samples can be specified via the group_by
    parameter.

    Groups appear as entries in the legend by default (i.e., when legend_factors
    is not True). What text is displayed for each group can be controlled using
    the make_group_label and labelers parameters. The labelers should be a list
    of functions of length no greater than group_by. labelers[i] should be a
    function that maps a value for column group_by[i] in the sample metadata to
    a value to use as the "label" for that value. If labelers[i] is None, or if
    i >= len(labelers), then the label for the value is just taken to be the
    value itself.

    make_group_label is a function that is used to construct the actual label
    that appears in the legend from the individual labels of the values in each
    of the columns for the group. The argument to the function is an iterable
    over the labels of the values of the columns for the group. The function
    should return a string to be used to represent the group in the legend. When
    legend_factors is not None, the function needs to be able to produce a label
    for only some of the columns; the iterable may have length less than that of
    group_by.

    draw_pcoa uses the colors parameter to select colors for the markers
    representing samples of a group. colors can either be a
    matplotlib.colors.ListedColormap, a sequence of RGB values, or a mapping
    from group values to RGB values. In the first two cases, groups are assigned
    colors in sorted order by their values, or, if order_by is set, in the order
    they appear in the sample_metadata DataFrame. In the last case, draw_pcoa
    looks up the group values in the mapping to determine the color.

    This function also provides a more sophisticated way of specifying how group
    values should be mapped to visual encodings. The index_to_kwargs,
    group_to_kwargs, and index_group_to_kwargs parameters allow the caller to
    provide arbitrary functions to map group values and indices to keyword
    arguments to be provided to ax.scatter.

    index_to_kwargs should map a group's index to a dict containing keyword
    arguments to pass to ax.scatter for drawing samples in that group. By
    default, a group's index is the group value's position when all group values
    are sorted. If order_by has been provided to draw_pcoa, then the group's
    index is the group value's position when all group values are sorted
    according to when they first appear in the sample_metadata DataFrame.

    group_to_kwargs should map a group's value to a dict containing keyword
    arguments to pass to ax.scatter for drawing the samples in that group.

    index_group_to_kwargs should map a group's index AND its value to a dict
    containing keyword arguments to pass to ax.scatter for the samples in that
    group.

    The order_by parameter can be used to sort the sample_metadata DataFrame by
    specified columns. As explained above, the order_by parameter changes how
    groups are assigned indices, indirectly affecting how groups are assigned
    colors and other visual encodings. Since groups appear in the legend based
    on their indices, order_by also affects the order in which groups appear in
    the legend.

    order_by can be used not only with the sample_metadata columns, but with the
    principal component coordinates. To sort by principal component i, provide
    f"PC{i}" as the order_by column.

    The sort_key parameter can be used to sort the sample_metadata DataFrame by
    values that do not appear explicitly in the columns of sample_metadata (or
    the principal component coordinates). sort_key should be a list of functions
    such that sort_key[i] maps column group_by[i] of the sample_metadata to a
    Pandas Series to use for sorting instead of column group_by[i].

    draw_pcoa can optionally draw an ellipsoid for each group of samples. This
    behavior can be enabled by setting the ellipsoids parameter to True. The
    specific ellipsoid to be drawn for each group can be controlled via the
    make_ellipsoid parameter. make_ellipsoid should be a function accepting the
    coordinates of all samples in a group and should return an
    confidence_ellipsoid.Ellipsoid object. Parameters for how the ellipsoids
    should be drawn can be specified via the ellipsoid_kwargs argument.

    By default, the relative contribution of a principal component axis is shown
    on the axis's label. The relative contribution of a component is defined as
    the eigenvalue of the component divided by the sum of the eigenvalues of all
    possible components (including those beyond the number of dimensions shown
    in the plot). To disable showing relative contribution on the axes, pass
    False for the contribution parameter.

    Individual samples can be annotated with their names (based on the
    sample_name_column) by passing True for the annotate parameter. By default,
    each sample label is placed directly on top of the scatter plot marker to
    which it refers. Since this behavior can sometimes cause labels to overlap
    and hurt readability, the adjust parameter is offered. When adjust is True,
    the adjustText library is used to automatically move labels to avoid
    overlaps. The points to which the labels refer is indicated via line
    segments between the points and the labels.

    By default, this function shows every group value as a separate entry in the
    legend, but it is sometimes desirable to "factor" the legend entries,
    showing how different variables affect different aspects of visual encoding
    independently. 

    When legend_factors is not None, for each sequence of columns specified in
    legend_factors, and each of the values for those columns, this function
    finds the scatter keyword argument key-value pairs shared among all groups
    that have those values in those columns. These common keyword arguments are
    assumed to be associated with groups having those specific values in those
    columns, and a legend entry is created to reflect that. For each sequence of
    columns in the legend factors, and each value for those columns for which
    some common scatter keyword argument key-value pairs were found, a legend
    entry is created to show that having those values for those columns produces
    a point marker created with all of those common key-value pairs as keyword
    arguments.

    Every marker must have certain visual properties, but not every set of
    common key-value pairs found by draw_pcoa will specify all of those
    properties. Any missing properties will be provided using the given
    default_legend_marker argument, which should be a dict containing default
    keyword arguments to ax.scatter.    

    Parameters:
        dis_df:                   DataFrame containing genetic dissimilarities.
        sample_metadata:          DataFrame containing sample metadata.
        group_by:                 Column to use as category for scatter plot.
        sample_name_column (str): Column to use as the sample name.
        make_group_label:         Function to make label from group name.
        labelers:                 Functions to apply to each group name.
        colors:                   "Mapping" from group name or index to color.
        legend (bool):            Whether to draw legend.
        index_to_kwargs:          Mapping from group index to scatter kwargs.
        group_to_kwargs:          Mapppng from group name to scatter kwargs.
        index_group_to_kwargs:    Mapping from name + index to scatter kwargs.
        order_by:                 Columns on which to order samples for legend.
        sort_key:                 Function to get keys for sorting samples.
        ellipsoids:               Draw ellipsoids on the plot.
        make_ellipsoid:           Function to obtain stat ellipsoid from data.
        ellipsoid_kwargs (dict):  kwargs to pass to draw_ellipsoid.
        contribution (bool):      Show relative contribition on axes.
        annotate (bool):          Label individual points on the plot.
        adjust (bool):            Try adjusting annotations for readability.
        legend_factors:           Separate independent encodings in legend.
        default_legend_marker:    Default kwargs for marker to use in legend.
        dropna (bool):            Drop rows with group keys having NA values.
        dimensions (int):         Number of dimensions (2 or 3) for the plot.
        ax:                       The matplotlib Axes on which to draw.

    """
    def make_label(x, l):
        x = list(x)
        for i in range(len(x)):
            try:
                x[i] = l[i](x[i])
            except TypeError:
                pass
        return make_group_label(x)
    if isinstance(group_by, str):
        group_by = [group_by]
    if ellipsoid_kwargs is None:
        ellipsoid_kwargs = {}
    if ellipsoids is True:
        ellipsoids = list(group_by)
    if labelers is None:
        labelers = [None]*len(group_by)
    if legend_factors is True:
        legend_factors = group_by
    if axis_label_kwargs is None:
        axis_label_kwargs = {}
    proj = dim_to_proj[dimensions]
    if ax is None:
        ax = plt.gca()
        if ax.name != proj:
            fig = ax.figure
            ax.remove()
            ax = fig.add_subplot(projection=proj)
    group_to_index = dict(map(reversed, enumerate(group_by)))
    dm = skb.DistanceMatrix(dis_df, ids=dis_df.columns)
    pcoa_results = skb.stats.ordination.pcoa(dm)
    pcs = [f"PC{x}" for x in range(1, dimensions + 1)]
    joined = sample_metadata.join(
        pcoa_results.samples[pcs],
        sample_name_column
    )
    # print(joined)
    ellipsoid_group_to_color = {}
    if order_by:
        if not isinstance(order_by, list):
            order_by = [order_by]
        if callable(sort_key):
            sort_key = [sort_key]
        if sort_key is not None:
            sort_key = [s if s is not None else id_ for s in sort_key]
        joined = _keyed_multi_sort(joined, order_by, sort_key)
    texts = []
    factored = defaultdict(dict)
    for i, (ge, df) in enumerate(
            joined.groupby(
                group_by,
                sort=not order_by,
                dropna=dropna
            )
    ):
        kwargs = dict(dimension_default_kwargs[dimensions])
        kwargs |= index_group_to_kwargs(i, ge)
        kwargs |= index_to_kwargs(i)
        try:
            kwargs |= group_to_kwargs(ge)
        except TypeError:
            kwargs |= group_to_kwargs(*ge)
        kwargs.setdefault("label", make_label(ge, labelers))
        color = None
        try:
            color = colors.colors[i]
        except AttributeError:
            try:
                color = colors[ge]
            except IndexError:
                colors = colors[i]
            except TypeError:
                pass
        if color is not None:
            kwargs.setdefault("color", color)
        kwargs |= scatter_kwargs
        #print(df, kwargs)
        # if kwargs["marker"] == "o":
        #     print("Scatter", [df[c] for c in pcs])
        ax.scatter(
            *(df[c] for c in pcs),
            **kwargs
        )
        if ellipsoids:
            new_kwargs = dict(ellipsoid_kwargs)
            try:
                new_kwargs.setdefault("color", kwargs["color"])
            except KeyError:
                pass
            key = tuple(ge[group_to_index[g]] for g in ellipsoids)
            try:
                if ellipsoid_group_to_color[key] != new_kwargs:
                    raise ValueError(
                        "All markers for an ellipsoid group must have the same "
                        "color."
                    )
            except KeyError:
                pass
            ellipsoid_group_to_color[key] = new_kwargs
        if annotate:
            for ix, r in df.set_index(sample_name_column).iterrows():
                texts.append(
                   ax.text(*(r[c] for c in pcs), ix, size="xx-small")
                )
        if legend_factors is not None:
            for t in map(as_tuple, legend_factors):
                key = tuple(ge[group_to_index[x]] for x in t)
                if key not in factored[t]:
                    factored[t][key] = set(kwargs.items())
                else:
                    factored[t][key] &= set(kwargs.items())
    if ellipsoids:
        for i, (ge, df) in enumerate(joined.groupby(ellipsoids, dropna=dropna)):
            if df.shape[0] > dimensions:
                draw_ellipsoid(
                    make_ellipsoid(df[pcs]),
                    ax=ax,
                    **ellipsoid_group_to_color[ge]
                )
    labels = list(pcs)
    if contribution:
        total_eig = sum(pcoa_results.eigvals)
        for i in range(len(labels)):
            labels[i] = "{} ({:05.2f}%)".format(
                labels[i],
                100 * pcoa_results.eigvals.iloc[i] / total_eig
            )
    for letter, la in zip("xyz", labels):
        getattr(ax, f"set_{letter}label")(la, **axis_label_kwargs)
    if annotate and adjust:
        adjust_text(
            texts,
            arrowprops=dict(arrowstyle="->", color='black', lw=0.5),
            # force_text = (0.3, 0.3)
        )
        pass
    if legend:
        if legend_factors is not None:
            handles = []
            names = []
            for cat, entries in factored.items():
                labeler = functools.partial(
                    make_label,
                    l=[labelers[group_to_index[x]] for x in cat]
                )
                for entry, kwargs in entries.items():
                    kw = default_legend_marker | dict(kwargs)
                    kw = {
                        scatter_to_line2d.get(k, k): v for (k, v) in kw.items()
                        if k not in scatter_ignore
                    }
                    handles.append(
                        mpl.lines.Line2D(
                            [0],
                            [0],
                            **kw
                        )
                    )
                    names.append(labeler(entry))
            ax.legend(handles, names)
        else:
            ax.legend()
    return pcoa_results

def dump_emperor_data(pcoa_results, sample_metadata, path):
    emp = Emperor(pcoa_results, sample_metadata, remote=".")
    with open(path, "w") as data_file:
        json.dump(
            emp._to_dict(
                emp._process_data(
                    emp.custom_axes,
                    emp.jackknifing_method)
            ),
            data_file,
            indent=2
        )

draw_pcoa_2d = functools.partial(draw_pcoa, dimensions=2)
draw_pcoa_3d = functools.partial(draw_pcoa, dimensions=3)
