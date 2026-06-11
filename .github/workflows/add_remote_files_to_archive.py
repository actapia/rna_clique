import argparse
import zipfile
import tarfile
import lzma
import magic
import gzip
import contextlib

import requests

from io import BytesIO
from pathlib import Path
from urllib.parse import urlparse


def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("--remote", type=urlparse, nargs="*", default=[])
    parser.add_argument("--merge-zips", "--merge-archives", action="store_true")
    parser.add_argument("--remote-dir")
    parser.add_argument("--compress-level", type=int, default=1)
    return parser.parse_args()

def get_files_zip(zip_file):
    with zipfile.ZipFile(zip_file, "r") as merge_file:
        for f in merge_file.infolist():
            if not f.is_dir():
                yield f.filename, merge_file.read(f)

def get_files_tar(tar_file):
    with tarfile.open(fileobj=tar_file, mode="r") as merge_file:
        for f in merge_file.getmembers():
            if f.isreg():
                yield f.path, merge_file.extractfile(f).read()

openers = {
    "application/zip": get_files_zip,
    "application/x-tar": get_files_tar
}

decompressors = {
    "application/gzip": lambda x, *a, **kw: gzip.GzipFile(*a, fileobj=x, **kw),
    "application/x-xz": lambda x, *a, **kw: lzma.LZMAFile(x, *a, **kw),
}

def main():
    args = handle_arguments()
    if args.remote_dir:
        remote_prefix = f"{args.remote_dir}"
    else:
        remote_prefix = ""
    with zipfile.ZipFile(
            args.archive,
            "a",
            compression=zipfile.ZIP_DEFLATED
    ) as archive:
        for url in args.remote:
            resp = requests.get(url.geturl())
            resp.raise_for_status()
            if args.merge_zips:
                bio = BytesIO(resp.content)
                with contextlib.ExitStack() as stack:
                    mime = resp.headers.get(
                        "Content-Type"
                    )
                    try:
                        bio = stack.enter_context(
                            decompressors[mime](bio, mode="rb")
                        )
                        mime = magic.from_buffer(bio.read(1024), mime=True)
                        bio.seek(0)                        
                    except KeyError:
                        pass
                    for path, contents in openers[mime](bio):
                        archive.writestr(
                            "{}/{}".format(
                                remote_prefix,
                                path
                            ),
                            contents,
                            compresslevel=args.compress_level
                        )
            else:
                archive.writestr(
                    "{}/{}".format(
                        remote_prefix,
                        Path(url.path).name
                    ),
                    resp.content,
                    compresslevel=args.compress_level,
                )
            
if __name__ == "__main__":
    main()

                
