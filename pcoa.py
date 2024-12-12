import skbio as skb
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt

from plots import default_group_label_maker, _keyed_multi_sort

from typing import Callable, Union, Optional
from collections.abc import Iterable, Sequence, Mapping

def empty_dict(*args, **kwargs):
    return {}

def draw_pcoa_2d(
        dis_df: pd.DataFrame,
        sample_metadata: pd.DataFrame,
        group_by: str | Iterable[str],
        sample_name_column: str = "name",
        make_group_label: Callable[[Iterable], str] = default_group_label_maker,
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
        **scatter_kwargs
) -> skb.stats.ordination.OrdinationResults:
    plt.rcParams["figure.figsize"] = (6.4, 4.8)
    dm = skb.DistanceMatrix(dis_df, ids=dis_df.columns)
    pcoa_results_2d = skb.stats.ordination.pcoa(dm, number_of_dimensions=2)
    joined = sample_metadata.join(
        pcoa_results_2d.samples[["PC1","PC2"]],
        sample_name_column
    )
    if order_by:
        if not isinstance(order_by, list):
            order_by = [order_by]
        if callable(sort_key):
            sort_key = [sort_key]
        joined = _keyed_multi_sort(joined, order_by, sort_key)
    for i, (ge, df) in enumerate(joined.groupby(group_by, sort=not order_by)):
        kwargs = index_group_to_kwargs(i, ge)
        kwargs |= index_to_kwargs(i)
        try:
            kwargs |= group_to_kwargs(ge)
        except TypeError:
            kwargs |= group_to_kwargs(*ge)
        kwargs.setdefault("label", make_group_label(ge))
        color = None
        try:
            color = colors.colors[i]
        except AttributeError:
            try:
                color = colors[ge]
            except IndexError:
                color = colors[i]
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
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    if legend:
        plt.legend()
    return pcoa_results_2d
