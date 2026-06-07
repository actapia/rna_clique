import pandas as pd

from itertools import zip_longest
from typing import Callable, TypeVar, Optional
from collections.abc import Iterable

T = TypeVar("T")

def default_group_label_maker(group_values: str | Iterable[str]) -> str:
    """Return a single str or return multiple values joined by commas."""
    if isinstance(group_values, str):
        return group_values
    else:
        return ", ".join(map(str, group_values))

def _keyed_multi_sort(
        df: pd.DataFrame,
        columns: list[str],
        keys: Optional[Iterable[Callable] | Callable] = None
) -> pd.DataFrame:
    """Sort a dataframe on multiple columns, optionally using keys for each one.

    This function sorts a Pandas dataframe by multiple columns. The
    DataFrame.sort_values method can do the same thing but does not support
    using different key functions for each column. This function adds that
    feature.

    This function should be provided a DataFrame and a list of columns. If
    provided, the keys parameter should be an Iterable (typically a list) of
    Callables or a single Callable. If the keys parameter is a list of
    Callables, it should have length no greater than the list of columns. The
    ith key in the keys Iterable is assumed to correspond to the ith elememnt of
    the columns list; the ith key will be used to sort the ith column. If no key
    should be used for a specific column, the correwsponding value in the keys
    Iterable can be None. If the keys Iterable is shorter than the columns
    list, then all columns not associated with corresponding key functions will
    be sorted without a key.

    If there is only one element in both the keys and columns, the function will
    simply call the sort_values method for the given column and key. Otherwise,
    this function will sort on each of the columns and corresponding keys in
    sequence using a stable sort.

    When a Callable is provided for the keys parameter instead of an Itearble of
    Callables, this function behaves like DataFrame.sort_values, using the same
    key for every column.

    Parameters:
        df:             DataFrame to sort.
        columns (list): Columns on which to sort.
        keys:           Keys for transforming values on which to sort.

    Returns:
        A dataframe with the same rows, sorted by the given columns and keys.
    """    
    if keys is None:
        return df.sort_values(columns)
    else:
        try:
            if len(columns) == 1 and len(keys) == 1:
                return df.sort_values(columns, key=keys[0])
        except TypeError:
            return df.sort_values(columns, key=keys)
        for column, key in reversed(list(zip_longest(columns, keys))):
            df = df.sort_values(column, key=key, kind="stable")
        return df

def _composite_transform(name: str) -> Callable[[T], T]:
    """Get a method that applies a composite transform.

    This function gets a method that applies a composite transform by applying
    each of the component transforms in sequence.

    Parameters:
        name (str): Name of transformation method to call on components.

    Returns:
        A composite transform method calling components' methods in sequence.
    """
    def inner(self, o):
        for t in self.transforms:
            o = getattr(t, name)(o)
        return o
    return inner

class BasicCompositeTransform:
    """A composition of multiple matplotlib transforms.
    
    This function applies multiple matplotlib transforms in sequence, composing
    them. Hence, the range of each function should be a subset of the domain of
    the next.

    Attributes:
        transforms: The transforms to apply, in order.
    """
    def __init__(self, *args):
        """Construct a BasicCompositeTransform with the given transforms.

        The arguments provided will be the transforms, in order, to apply.
        """
        self.transforms = args

    def inverted(self):
        """Return the inverse transform for this transform.

        The inverse of the composite function is the transform containing the
        inverses of each of the components, in reverse order.

        The inverse should have the property that, for any value x in the domain
        of a transform t,
        BasicCompositeTransform(t, t.inverted()).transform_point(t) == t.

        Returns:
            The inverse transform of this transform.
        """
        return BasicCompositeTransform(
            *(
                t.inverted()
                for t in reversed(self.transforms)
            )
        )
            
        
    # def transform_point(self, p):
    #     """Transform the point by applying each transform in order."""
    #     for t in self.transforms:
    #         p = t.transform_point(p)
    #     return p

    # def transform_bbox(self, p):
    #     """Transform the Bbox by applying each transform in order."""

# This automatically creates the transform_point and transform_bbox methods of
# BasicCompositeTransform.
for n in ["point", "bbox"]:
    n = f"transform_{n}"
    setattr(BasicCompositeTransform, n, _composite_transform(n))

def _transform_ax(trans, coords: Iterable[float], axis: int):
    """Transform points on a single axis.

    Parameters:
        trans:      Transformation that can transform multidimensional points.
        coords:     Iterable of points on a single axis to transform.
        axis (int): Axis on which to perform the transformation.

    Returns:
        The points, transformed on the specified axis.
    """
    def get_point(coord: float) -> tuple[float, float]:
        p = [0, 0]
        p[axis] = coord
        return p
    return [trans.transform_point(get_point(coord))[axis] for coord in coords]

def as_tuple(x) -> tuple:
    """Return x as a tuple, converting or creating a singleton."""
    if isinstance(x, str):
        return (x,)
    else:
        try:
            return tuple(x)
        except TypeError:
            return (x,)
