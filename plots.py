from itertools import zip_longest

from collections.abc import Iterable

def default_group_label_maker(group_values: str | Iterable[str]) -> str:
    """Return a single str or return multiple values joined by commas."""
    if isinstance(group_values, str):
        return group_values
    else:
        return ", ".join(map(str, group_values))

def _keyed_multi_sort(df, columns, keys=None):
    if keys is None:
        return df.sort_values(columns)
    elif len(columns) == 1 and len(keys) == 1:
        return df.sort_values(columns, key=keys[0])
    for column, key in reversed(list(zip_longest(columns, keys))):
        df = df.sort_values(column, key=key, kind="stable")
    return df

class BasicCompositeTransform:
    """A composition of multiple matplotlib transforms."""
    def __init__(self, *args):
        self.transforms = args
        
    def transform_point(self, p):
        """Transform the point by applying each transform in order."""
        for t in self.transforms:
            p = t.transform_point(p)
        return p

def _transform_ax(trans, coords, axis):
    def get_point(coord):
        p = [0, 0]
        p[axis] = coord
        return p
    return [trans.transform_point(get_point(coord))[axis] for coord in coords]
