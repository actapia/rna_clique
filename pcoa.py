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

from confidence_ellipsoid import Ellipsoid, draw_ellipse, conf_ellipsoid
from plots import as_tuple
from similarity_computer import id_

def empty_dict(*args, **kwargs):
    return {}

default_marker_style = {"color": "black", "ls": "", "marker": "d"}
scatter_to_line2d = {
    "edgecolors": "mec",
    "linewidth": "mew",
    "s": "ms"
}

def draw_pcoa_2d(
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
        ellipses: bool | Iterable[str] = False,
        make_ellipse: Callable[
            [pd.DataFrame],
            Ellipsoid
        ] = conf_ellipsoid(0.95),
        ellipse_kwargs: dict[str, Any] = None,
        contribution: bool = True,
        annotate: bool = False,
        legend_factors: Optional[list[tuple] | bool] = None,
        default_legend_marker: dict = default_marker_style,
        #coordinate_sort: bool = False,
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
        ellipses:                 Draw ellipses on the plot.
        make_ellipse:             Function to obtain stat ellipse from data.
        ellipse_kwargs (dict):    kwargs to pass to draw_ellipse.
        contribution (bool):      Show relative contribition on axes.
        annotate (bool):          Label individual points on the plot.
        legend_factors:           Separate independent encodings in legend.
        default_legend_marker:    Default kwargs for marker to use in legend.
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
    if ellipse_kwargs is None:
        ellipse_kwargs = {}
    if ellipses is True:
        ellipses = list(group_by)
    if labelers is None:
        labelers = [None]*len(group_by)
    if legend_factors is True:
        legend_factors = group_by
    group_to_index = dict(map(reversed, enumerate(group_by)))
    dm = skb.DistanceMatrix(dis_df, ids=dis_df.columns)
    pcoa_results = skb.stats.ordination.pcoa(dm)
    joined = sample_metadata.join(
        pcoa_results.samples[["PC1","PC2"]],
        sample_name_column
    )
    ellipse_group_to_color = {}
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
    for i, (ge, df) in enumerate(joined.groupby(group_by, sort=not order_by)):
        kwargs = index_group_to_kwargs(i, ge)
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
        plt.scatter(
            df["PC1"],
            df["PC2"],
            **kwargs
        )
        if ellipses:
            new_kwargs = dict(ellipse_kwargs)
            try:
                new_kwargs.setdefault("color", kwargs["color"])
            except KeyError:
                pass
            key = tuple(ge[group_to_index[g]] for g in ellipses)
            try:
                if ellipse_group_to_color[key] != new_kwargs:
                    raise ValueError(
                        "All markers for an ellipse group must have the same "
                        "color."
                    )
            except KeyError:
                pass
            ellipse_group_to_color[key] = new_kwargs
        if annotate:
            for ix, r in df.set_index("run_accession").iterrows():
                texts.append(
                    plt.gca().text(r["PC1"], r["PC2"], ix, size="xx-small")
                )
        if legend_factors is not None:
            for t in map(as_tuple, legend_factors):
                key = tuple(ge[group_to_index[x]] for x in t)
                if key not in factored[t]:
                    factored[t][key] = set(kwargs.items())
                else:
                    factored[t][key] &= set(kwargs.items())
    if ellipses:
        for i, (ge, df) in enumerate(joined.groupby(ellipses)):
            draw_ellipse(
                make_ellipse(df[["PC1", "PC2"]]),
                **ellipse_group_to_color[ge]
            )
    labels = ["PC1", "PC2"]
    if contribution:
        total_eig = sum(pcoa_results.eigvals)
        for i in range(len(labels)):
            labels[i] = "{} ({:05.2f}%)".format(
                labels[i],
                100 * pcoa_results.eigvals.iloc[i] / total_eig
            )
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
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
            plt.legend(handles, names)
        else:
            plt.legend()
    return pcoa_results
