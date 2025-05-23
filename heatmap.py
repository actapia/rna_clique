import functools
import typing

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.transforms
from matplotlib import pyplot as plt

from plots import (
    default_group_label_maker,
    _keyed_multi_sort,
    BasicCompositeTransform,
    _transform_ax
)

from collections.abc import Sequence
from typing import Union, Optional, Callable, Literal, Any
from itertools import zip_longest
from collections.abc import Iterable

#from IPython import embed

# noinspection PyTypeChecker
axis_to_pos = dict(map(reversed, enumerate(["x", "y"])))
axis_to_ha = {"y": "left", "x": "center"}
axis_to_va = {"y": "center", "x": "bottom"}

def draw_heatmap(
        mat: pd.DataFrame,
        *heatmap_args,
        sample_metadata: Optional[pd.DataFrame] = None,
        sample_name_column: str = "name",
        order_by: Optional[Union[str, Iterable[str]]] = None,
        square: bool = True,
        group_by: Optional[Union[str, Iterable[str]]] = None,
        draw_group_labels: bool = False,
        make_group_label: Callable[[Iterable], str] = default_group_label_maker,
        digit_annot: Optional[int] = None,
        label_padding_x: float = 0.0275,
        label_padding_y: float = 0.0275,
        sort_key: Optional[
            Union[
                Callable[[pd.Series], pd.Series],
                Sequence[Callable[[pd.Series], pd.Series]]
            ]
        ] = None,
        label_kwargs: Optional[dict[str, Any]] = None,
        x_label_kwargs: Optional[dict[str, Any]] = None,
        y_label_kwargs: Optional[dict[str, Any]] = None,
        draw_debug_points: bool = False,
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
        group_by:                 Columns to use for grouping.
        draw_group_labels (bool): Whether to draw labels for order_by values.
        make_group_label:         Function to get group label from group values.
        digit_annot (int):        Annotate with the most significant digits.
        label_padding_x (float):  x padding to add to group labels on y axis.
        label_padding_y (float):  y padding to add to group labels on x axis.
        sort_key:                 Function to transform sort column to sort key.
        label_kwargs (dict):      kwargs to give to plt.text for group labels.
        x_label_kwargs (dict):    kwargs to give to plt.text for x group labels.
        y_label_kwargs (dict):    kwargs to give to plt.text for y group labels.
        draw_debug_points (bool): Whether to plot points useful for debugging.
    """

    def _draw_group_labels(
            axis: Literal["x", "y"],
            padding: float = 0,
            **kwargs
    ):
        
        def draw_labels(ppos=0):
            texts = []
            max_width = 0
            for group, df in sample_metadata.groupby(group_by):
                avg_pos = sum(
                    pos_dict[sample] for sample in df["index"]
                )/len(df)
                coords = [0, 0]
                coords[para] = avg_pos
                coords[perp] = ppos
                group_label = make_group_label(typing.cast(Iterable, group))
                # PyCharm's type checker doesn't seem to recognize that *pos
                # covers the first two arguments.
                # noinspection PyTypeChecker
                texts.append(
                    plt.text(
                        *coords,
                        group_label,
                        horizontalalignment=axis_to_ha[axis],
                        verticalalignment=axis_to_va[axis],
                        **kwargs
                    )
                )
                # canvas.renderer exists, but only when we're actually drawing.
                # noinspection PyUnresolvedReferences
                max_width = max(
                    max_width,
                    abs(
                        texts[-1].get_window_extent(
                            plt.gcf().canvas.renderer
                        ).transformed(ax.transData.inverted()).bounds[2+perp]
                    )
                )
            return texts, max_width

        tick_labels = getattr(ax, f"get_{axis}ticklabels")()
        para = axis_to_pos[axis]
        perp = (para + 1) % 2
        pos_dict = {
            t.get_text(): t.get_position()[para]
            for t in tick_labels
        }
        sign = para*2-1
        # noinspection PyUnresolvedReferences
        edge = sign*min(
            min(
                sign*t.get_window_extent(
                    plt.gcf().canvas.renderer
                ).transformed(ax.transData.inverted()).get_points()[:, perp]
            )
            for t in tick_labels
        )
        debug_point = [0, 0]
        debug_point[perp] = edge
        for t in pos_dict.values():
            debug_point[para] = t
            if draw_debug_points:
                plt.scatter(*debug_point, clip_on=False, color="red")
        if draw_debug_points:
            for t in tick_labels:
                for point in t.get_window_extent(
                        plt.gcf().canvas.renderer
                ).transformed(ax.transData.inverted()).get_points():
                    plt.scatter(*point, clip_on=False, color="blue", s=10)
        text_elems, m_width = draw_labels()
        for text in text_elems:
            text.remove()
        a2d = ax_to_data([padding], perp)[0] - ax_to_data([0], perp)[0]
        perp_pos = edge - sign*m_width - a2d
        draw_labels(perp_pos)
        return perp_pos
    
    if order_by:
        if isinstance(order_by, str):
            order_by = [order_by]
        sample_metadata = _keyed_multi_sort(
            sample_metadata.set_index(
                sample_name_column
            ).loc[mat.columns].reset_index(),
            order_by,
            keys=sort_key
        )
        #embed()
        mat = mat.iloc[sample_metadata.index].iloc[:, sample_metadata.index]
    if group_by:
        if isinstance(group_by, str):
            group_by = [group_by]
        assert order_by[:len(group_by)] == group_by
    elif order_by:
        group_by = order_by
    if digit_annot is not None:
        heatmap_kwargs["annot"] = (
            mat * 10**((-np.floor(np.log10(mat))).min(None) + digit_annot - 1)
        ).round()
    if label_kwargs is None:
        label_kwargs = {}
    if x_label_kwargs is None:
        x_label_kwargs = {}
    if y_label_kwargs is None:
        y_label_kwargs = {}
    label_kwargs.setdefault("fontsize", 10)
    x_label_kwargs = label_kwargs | x_label_kwargs
    y_label_kwargs = label_kwargs | y_label_kwargs
    # This trick comes from PyRsquared on Stack Overflow
    # https://stackoverflow.com/a/49608671
    mask = np.zeros_like(mat)
    mask[np.tril_indices_from(mask)] = True
    sns.set(rc={"figure.figsize": (15, 7)})
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
    ax.apply_aspect()
    fig = ax.get_figure()
    plt.sca(ax)
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
    major_kwargs = base_kwargs | {"lw": major_weight}
    minor_kwargs = base_kwargs | {"lw": minor_weight}
    fig.canvas.draw()
    ax_to_data = functools.partial(
        _transform_ax,
        BasicCompositeTransform(ax.transAxes, ax.transData.inverted())
    )
    data_to_inches_trans = BasicCompositeTransform(
            ax.transData,
            fig.dpi_scale_trans.inverted()
    )
    data_to_inches = functools.partial(
        _transform_ax,
        data_to_inches_trans,
    )
    inches_to_figure_trans = BasicCompositeTransform(
        fig.dpi_scale_trans,
        fig.transFigure.inverted()
    )
    inches_to_figure = functools.partial(
        _transform_ax,
        inches_to_figure_trans
    )
    figure_to_inches = functools.partial(
        _transform_ax,
        inches_to_figure_trans.inverted()
    )
    group_label_dist_x = 0
    group_label_dist_y = ax_to_data([1], 1)[0]
    if draw_group_labels:
        group_label_dist_x = _draw_group_labels(
            "y",
            label_padding_x,
            **y_label_kwargs
        )
        group_label_dist_y = _draw_group_labels(
            "x",
            label_padding_y,
            **x_label_kwargs
        )
    xpos = [group_label_dist_x] + ax_to_data([1], 0)
    ypos = [0] + [group_label_dist_y]
    if group_by:
        for pos in sample_metadata.reset_index(
                drop=True
        ).groupby(group_by[0]).tail(1).head(-1).index:
            pos += 1
            ax.plot(xpos, [pos, pos], clip_on=False, **major_kwargs)
            ax.plot([pos, pos], ypos, clip_on=False, **major_kwargs)
            if draw_debug_points:
                ax.scatter([pos], [ypos[0]], clip_on=False, color="yellow")
                ax.scatter([pos], [ypos[1]], clip_on=False, color="purple")


        if len(group_by) > 1:
            for pos in sample_metadata.reset_index(
                    drop=True
            ).groupby(order_by).tail(1).head(-1).index:
                pos += 1
                ax.plot(xpos, [pos, pos], clip_on=False, **minor_kwargs)
                ax.plot([pos, pos], ypos, clip_on=False, **minor_kwargs)
        # noinspection PyTypeChecker
        plt.xlabel(None)
        # noinspection PyTypeChecker
        plt.ylabel(None)
        new_mins = [
            min(
                data_to_inches([group_label_dist_x], 0)[0],
                fig.bbox_inches.extents[0]
            ),
            min(
                data_to_inches([group_label_dist_y], 1)[0],
                fig.bbox_inches.extents[1]
            )
        ]
        shift = [0, 0]
        for i in range(len(new_mins)):
            if new_mins[i] < 0:
                shift[i] = inches_to_figure([-new_mins[i]], i)[0]
                new_mins[i] = 0
        for x in fig.axes:
            fig.canvas.draw()
            current_position = x.get_position(original=True)
            new_position = matplotlib.transforms.Bbox(
                    [
                        [
                            current_position.extents[0] + shift[0],
                            current_position.extents[1] + shift[1]
                        ],
                        [
                            current_position.extents[2] + shift[0],
                            current_position.extents[3] + shift[1]
                        ]
                    ]
                )
            x._set_position(new_position)
        fig.bbox_inches = matplotlib.transforms.Bbox(
            [
                new_mins,
                [
                    fig.bbox_inches.extents[2] + figure_to_inches(
                        [shift[0]], 0
                    )[0],
                    fig.bbox_inches.extents[3] + figure_to_inches(
                        [shift[1]], 1
                    )[0],
                ]
            ]
        )
    #fig.bbox = fig.bbox._transform.transform_bbox(fig.bbox_inches)
