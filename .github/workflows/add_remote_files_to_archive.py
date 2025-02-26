import argparse
import zipfile

import requests

from io import BytesIO
from pathlib import Path
from urllib.parse import urlparse


def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("--remote", type=urlparse, nargs="*", default=[])
    parser.add_argument("--merge-zips", action="store_true")
    parser.add_argument("--remote-dir")
    parser.add_argument("--compress-level", type=int, default=1)
    return parser.parse_args()

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
            if args.merge_zips and \
               resp.headers.get("Content-Type") == "application/zip":
                zip_bytes = BytesIO(resp.content)
                with zipfile.ZipFile(zip_bytes, "r") as merge_file:
                    for f in merge_file.infolist():
                        if not f.is_dir():
                            archive.writestr(
                                "{}/{}".format(
                                    remote_prefix,
                                    f.filename
                                ),
                                merge_file.read(f),
                                compresslevel=args.compress_level,
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

                
