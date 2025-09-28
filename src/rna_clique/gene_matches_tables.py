import itertools

import pandas as pd

from pathlib import Path
from collections.abc import Iterable
from typing import Iterator

def read_table(
        path: Path,
        head: int = None,
        head_unsupported: bool = True
) -> pd.DataFrame:
    """Read a dataframe from a specified path, guessing format from extension.

    The optional head parameter may be used to obtain only the first head rows
    of the table. For tables stored in the pickle format, this is accomplished
    by reading the whole table and taking the head of the table, so it's not any
    faster than reading the whole file. For tables stored in the HDF5 format, it
    is possible to read only the first head rows, which can be much faster than 
    reading the full file.

    Parameters:
        path:                    Path to file from which to read dataframe.
        head (int):              If specified, only head rows will be provided.
        head_unsupported (bool): Head dataframe even if unsupported by format.

    Returns:
        The dataframe from the file at the specified path.
    """
    if path.suffix == ".pkl":
        res = pd.read_pickle(path)
    elif path.suffix == ".h5":
        res = pd.read_hdf(path, stop=head)
    else:
        raise ValueError(
            f"Could not determine file type for extension {path.suffix}."
        )
    if head and head_unsupported:
        res = res.head(head)
    return res

def write_table(df: pd.DataFrame, path: Path):
    """Save a dataframe to a specified path, guessing format based on extension.

    Parameters:
        df:   The dataframe to save.
        path: Path to which data will be saved.
    """
    if path.suffix == ".pkl":
        df.to_pickle(path)
    elif path.suffix == ".h5":
        df.to_hdf(path, key="gene_matches", format="table")
    else:
        raise ValueError(
            f"Could not determine file type for extension {path.suffix}."
        )

def multi_glob(path: Path, globs: Iterable[str]) -> Iterator[Path]:
    """Return an iterator over multiple glob patterns.

    Parameters:
        path:  The path at which to glob.
        globs: The glob patterns to use.
    """
    return itertools.chain(*map(path.glob, globs))

def get_table_files(path: Path) -> Iterator[Path]:
    """Get stored gene matches talbes from a directory.

    This function simply uses glob patterns to look for files that could be
    serialized gene matches tables---it makes no effort to verify the contents
    of the files. If your directory contains non-gene matches table pickle or
    HDF5 files, this function will unwittingly get those files, too.

    Untrusted files should not be loaded this way since pickles allow arbitrary
    code exection. This may be true even for HDF5 files, which can indirectly
    cause pickles to be loaded.

    Parameters:
        path: Path to the directory containing the gene matches tables.
    """    
    return multi_glob(path, ["*.pkl", "*.h5"])
