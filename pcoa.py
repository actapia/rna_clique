import json
import functools
import skbio as skb
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
from adjustText import adjust_text

from plots import default_group_label_maker, _keyed_multi_sort

from collections import defaultdict
from typing import Callable, Union, Optional, Any
from collections.abc import Iterable, Sequence, Mapping

from confidence_ellipsoid import Ellipsoid, draw_ellipsoid, conf_ellipsoid
from plots import as_tuple
from similarity_computer import id_

def empty_dict(*args, **kwargs):
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

dim_to_proj = {2: "rectilinear", 3: "3d"}
dimension_default_kwargs = {2: {}, 3: {"depthshade": False}}

def draw_pcoa(
        dis_df: pd.DataFrame,
        sample_metadata: pd.DataFrame,
        group_by: str | Iterable[str],
        sample_name_column: str = "name",
        make_group_label: Callable[[Iterable], str] = default_group_label_maker,
        labelers: Optional[list[Callable[[Iterable], str]]] = None,
        colors: Optional[
            Union[
                mpl.colors.ListedColormap,
                Sequence[tuple[int, int, int]],
                Mapping[str, tuple[int, int, int]]
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
        legend_factors: Optional[list[tuple] | bool] = None,
        default_legend_marker: dict = default_marker_style,
        dropna: bool = True,
        dimensions: int = 2,
        ax: Optional[mpl.axes.Axes] = None,
        #coordinate_sort: bool = False,
        axis_label_kwargs: Optional[dict] = None,
        **scatter_kwargs
) -> skb.stats.ordination.OrdinationResults:
    """Draw a 2D PCoA plot and return the PCoA results.

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
    if annotate:
        adjust_text(
            texts,
            arrowprops=dict(arrowstyle="->", color='black', lw=0.5),
            force_text = (0.3, 0.3)
        )
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
