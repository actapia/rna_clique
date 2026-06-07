import functools
import tempfile
import io
import sys

import pandas as pd


from pathlib import Path
from typing import Callable
from collections.abc import Iterable

from . import config as config_module
from .app import set_except_hook

def write_hdf(df: pd.DataFrame, f: io.BytesIO, key: str = "matrix"):
    """Write the dataframe in HD5 format to the file-like object.

    Parameters:
        df:        The dataframe to write.
        f:         The file-like object to which to write the dataframe.
        key (str): The key under which to store the dataframe.
    """
    with pd.HDFStore(
            tempfile.mktemp(),
            driver="H5FD_CORE",
            driver_core_backing_store=0
    ) as store:
        store.put("matrix", df, errors="string", encoding="UTF-8")
        f.write(store._handle.get_file_image())

def ignore_kwargs(f: Callable, ignore: Iterable[str]) -> Callable:
    """Wrap a function, filtering certain kwargs.

    Parameters:
        f:      function to wrap
        ignore: kwargs to ignore

    Returns:
        A wrapped version of f that filters the specified kwargs.
    """
    ignore = set(ignore)
    def inner(*args, **kwargs):
        kwargs = {k: v for (k, v) in kwargs.items() if k not in ignore}
        return f(*args, **kwargs)
    return inner

ignore_header = functools.partial(ignore_kwargs, ignore={"header"})

# Functions for writing matrices in different formats.
writers = {
    "matrix": ignore_header(
        functools.partial(
            pd.DataFrame.to_csv,
            sep=" ",
            header=False,
            index=False
        ),
    ),
    "table": functools.partial(pd.DataFrame.to_csv, sep=" "),
    "csv": pd.DataFrame.to_csv,
    "hdf": ignore_header(write_hdf),
    "pickle": ignore_header(pd.DataFrame.to_pickle),
}

extension_to_format = {
    "h5": "hdf",
    "pkl": "pickle",
    "csv": "csv"
}

def build_parser():
    arg_config = config_module.RNACliqueConfigArgumentManager(
        description="Export a computed dissimilarity matrix.",
    )
    arg_config.expose_fields_with_default_aliases(
        "matrix",        
        required=True
    )
    arg_config.expose_fields_with_default_aliases(
        "output_dir",
    )
    arg_config.add_argument(
        "--export-out",
        "-x",
        type=Path,
        help="Path to which to export the matrix."
    )
    arg_config.add_argument(
        "-f",
        "--format",
        choices=writers,
        default={
            ("export_out",): lambda export_out: extension_to_format.get(
                export_out.suffix[1:],
                "matrix"
            ),
            (): "matrix",
        },
        help="Format for writing distance matrix."
    )
    arg_config.add_argument(
        "--header",
        action="store_true",
        help="Include header in distance matrix written."
    )
    return arg_config

def main():
    with set_except_hook():
        _, args, config = build_parser().get_arguments_and_config()
    with set_except_hook(config.verbose):
        mat = pd.read_hdf(config.matrix)
        if args.export_out:
            with open(args.export_out, "wb") as out:
                writers[args.format](mat, out, header=args.header)
        else:
            writers[args.format](mat, sys.stdout.buffer, header=args.header)

if __name__ == "__main__":
    main()
