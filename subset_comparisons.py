import os
import re

from pathlib import Path

import pandas as pd

from typing import Optional, Callable, Iterator
from collections.abc import Container, Iterable

from gene_matches_tables import read_table

default_filter_regex = re.compile("(.*)")

def matcher(
        included: Optional[Container[str]] = None,
        excluded: Optional[Container[str]] = None,
        include_regex: Optional[re.Pattern] = None
) -> Callable[[str], bool]:
    """Returns a function that checks if a string meets certain criteria.

    Specifically, the returned function returns a bool indicating whether its
    argument is in the provided Container or matches the given regex.

    Parameters:
        included:      A container of strings to be included.
        excluded:      A container of strings to be excluded.
        include_regex: A regular expression to match for inclusion.

    Returns:
        A function that checks if a str matches the regex or included strings.
    """
    def inner(x):
        return (
            (included is None and include_regex is None) or \
            (included is not None and x in included) or \
            (include_regex is not None and bool(include_regex.search(x)))
        ) and (excluded is None or x not in excluded)
    return inner

def relative_to(p1: Path, p2: Path) -> Path:
    """Returns the first path relative to the second."""
    return Path(os.path.relpath(str(p1), str(p2)))

def make_subset_comparisons(
        inputs: Iterable[Path],
        output_dir: Path,
        matches: Callable[[str], bool],
        sample_name_regex: re.Pattern
) -> Iterator[pd.DataFrame]:
    """Creates symlinks to stored dataframes whose samples satisfy a predicate.

    Parameters:
        inputs:            The Paths to the input dataframe pickles.
        output_dir:        The directory in which to create the symlinks.
        matches:           A function indicating whether a sample is included.
        sample_name_regex: A regular expression for parsing sample names.

    Returns:
        A generator yielding the dataframes whose samples satisfy the predicate.
    """
    for df_path in inputs:
        df = read_table(df_path, head=1, head_unsupported=False)
        if all(
                matches(
                    sample_name_regex.search(
                        Path(df[x + "sample"].iloc[0]).name
                    ).group(1)
                )
                for x in ["q", "s"]
        ):
            # We only need to re-read if it looks like we headed the table the
            # first time.
            if df.shape[0] == 1:
                df = read_table(df_path)
            dest = output_dir / df_path.name
            dest.symlink_to(relative_to(df_path, dest.parent))
            yield df

def handle_filters(include: Iterable[str], include_file: Path) -> set[str]:
    """Creates a set including both the given iterable and file's contents."""
    include = set(include)
    try:
        with open(include_file, "r") as filter_file:
            include |= set(l.rstrip() for l in filter_file)
    except TypeError:
        pass
    return include
