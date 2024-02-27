import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

from typing import Union, Optional
from collections.abc import Iterable

from IPython import embed

def draw_heatmap(
        mat: pd.DataFrame,
        *heatmap_args,
        sample_metadata: Optional[pd.DataFrame] = None,
        sample_name_column: str = "name",
        order_by: Optional[Union[str, Iterable[str]]] = None,
        square: bool = True,
        **heatmap_kwargs
):
    """Draw a heatmap representing a (dis)similarity matrix.

    This function uses sns.heatmap but applies a few changes to the plot that
    are appropriate for visualizing distance matrices. Any additional positional
    or keyword arguments provided to this function will be passed to
    sns.heatmap.

    Parameters:
        mat:                      (Dis)similarity matrix to visualize.
        sample_metadata:          Pandas dataframe containing sample metadata.
        sample_name_column (str): Column in sample_metadata for sample name.
        order_by:                 Columns by which to order samples.
        square (bool):            Whether to draw heatmap cells as squares.
    """    
    if order_by:
        if isinstance(order_by, str):
            order_by = [order_by]
        sample_metadata = sample_metadata.set_index(
            sample_name_column
        ).loc[mat.columns].reset_index().sort_values(order_by)
        #embed()
        mat = mat.iloc[sample_metadata.index].iloc[:,sample_metadata.index]
    # This trick comes from PyRsquared on Stack Overflow
    # https://stackoverflow.com/a/49608671
    mask = np.zeros_like(mat)
    mask[np.tril_indices_from(mask)] = True
    sns.set(rc={"figure.figsize": (15,7)})
    sns.set_style("white")
    ax = sns.heatmap(
        mat,
        *heatmap_args,
        xticklabels=mat.columns,
        yticklabels=mat.index,
        mask=mask,
        square=square,
        **heatmap_kwargs
    )
    yticks = plt.yticks(
        ha="right",
        multialignment="center",
        rotation=0,
        fontsize=7
    )
    plt.yticks(yticks[0] + 0.5, minor=True)
    xticks = plt.xticks(fontsize=6)
    plt.xticks(xticks[0] + 0.5, minor=True)
    ax.set_axisbelow(False)
    ax.grid(which="minor", alpha=0.15)
    divider_color = "0.8"
    major_weight = 2
    minor_weight = 1.5
    base_kwargs = {"color": divider_color, "alpha": 0.15}
    major_kwargs = base_kwargs |  {"lw": major_weight}
    minor_kwargs = base_kwargs |  {"lw": minor_weight}
    for pos in sample_metadata.reset_index(
            drop=True
    ).groupby(order_by[0]).tail(1).index:
        pos += 1
        plt.axvline(x=pos, **major_kwargs)
        plt.axhline(y=pos, **major_kwargs)
    if len(order_by) > 1:
        for pos in sample_metadata.reset_index(
                drop=True
        ).groupby(order_by[:2]).tail(1).index:
            pos += 1
            plt.axvline(x=pos, **minor_kwargs)
            plt.axhline(y=pos, **minor_kwargs)
    plt.xlabel(None)
    plt.ylabel(None)
