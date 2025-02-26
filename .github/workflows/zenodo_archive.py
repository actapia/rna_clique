import argparse
import getpass
import functools
import sys
import re
import os
import json
import datetime

from jsonpath_ng import jsonpath, parse

import requests

from pathlib import Path
try:
    from decrypt_token import decrypt_token
    from IPython import embed
except ImportError:
    pass

eprint = functools.partial(print, file=sys.stderr)

def add_environ_arg(parser, env, *args, **kwargs):
    parser.add_argument(
        *args,
        default=os.environ.get(env),
        required=env in os.environ,
        **kwargs
    )

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="*", type=Path)
    parser.add_argument("--sandbox", "-s", action="store_true")
    parser.add_argument("--token-file", "-t", type=Path)
    parser.add_argument("--print-vars", nargs="?", const=True, choices=["exit"])
    parser.add_argument("--clean", "-c", action="store_true")
    parser.add_argument("--publish", "-p", action="store_true")
    return parser.parse_args()

sandbox_zenodo_endpoint = "https://sandbox.zenodo.org/api/"
zenodo_endpoint = "https://zenodo.org/api/"

    
class ZenodoManager:
    methods = ["get", "head", "post", "put", "delete", "options", "patch"]
    
    def __init__(self, token, sandbox=False):
        self.token = token
        if sandbox:
            self.api_root = sandbox_zenodo_endpoint
        else:
            self.api_root = zenodo_endpoint

    def send_request(self, req_type, endpoint, *args, **kwargs):
        req = getattr(requests, req_type)
        kwargs.setdefault("headers", {})
        kwargs["headers"]["Authorization"] = f"Bearer {self.token}"
        res = req(self.api_root + endpoint, *args, **kwargs)
        res.raise_for_status()
        return res

    def get_paginated(self, *args, page_size=25, json_path="$", **kwargs):
        pattern = parse(json_path)
        kwargs.setdefault("params", {})
        kwargs["params"]["size"] = page_size
        kwargs["params"]["page"] = 1
        resp = self.get(*args, **kwargs)
        result = pattern.find(resp.json())[0].value
        total = result
        while result:
            #print(result)
            kwargs["params"]["page"] += 1
            resp = self.get(*args, **kwargs)
            result = pattern.find(resp.json())[0].value
            total = total + result
        return total

    def list_depositions(self, **params):
        return self.get_paginated("deposit/depositions", params=params)

    def get_deposition(self, id_):
        return self.get(f"deposit/depositions/{id_}")

    def create_deposition(self, data=None):
        if data is None:
            data = {}
        return self.post("deposit/depositions", json=data)

    def update_deposition(self, id_, data):
        return self.put(f"deposit/depositions/{id_}", json=data)

    def new_deposition_version(self, id_):
        return self.post(f"deposit/depositions/{id_}/actions/newversion")

    def upload_data(self, bucket, name, fp):
        #return self.put(f"{bucket}/{name}", data=fp)
        headers = {"Authorization": f"Bearer {self.token}"}
        return requests.put(f"{bucket}/{name}", data=fp, headers=headers)

    def upload_file(self, bucket, path: Path):
        with open(path, "rb") as fp:
            return self.upload_data(bucket, str(path), fp)

    def list_files(self, id_):
        return self.get(f"deposit/depositions/{id_}/files")

    def delete_file(self, id_, file_id):
        return self.delete(f"deposit/depositions/{id_}/files/{file_id}")

    def publish_deposition(self, id_):
        return self.post(f"deposit/depositions/{id_}/actions/publish")

for req_type in ZenodoManager.methods:
    setattr(
        ZenodoManager,
        req_type,
        (
            lambda r:
            lambda self, *args, **kwargs: self.send_request(r, *args, **kwargs)
        )(req_type)
    )

elasticsearch_re = re.compile("([:/])")
elasticsearch_escape = functools.partial(elasticsearch_re.sub, r"\\\1")

def jsonpath_find_one_if_exists(expr, x):
    res = expr.find(x)
    if len(res) == 0:
        return None
    elif len(res) > 1:
        raise ValueError("More than one match.")
    else:
        return res[0].value

environ_vars = [
    "GITHUB_TAG_URL",
    "GITHUB_URL",
    "PUBLISHED_AT",
    "RELEASE_DESCRIPTION",
    "REPO_NAME",
    "VERSION",
    "ZENODO_JSON",
]

def main():
    args = handle_arguments()
    if args.print_vars:
        for v in environ_vars:
            print("{}: {}".format(v, os.environ.get(v)))
    if args.print_vars == "exit":
        sys.exit(0)
    if args.token_file:
        with open(args.token_file, "rb") as token_file:
            token = decrypt_token(getpass.getpass(), token_file.read().rstrip())
    else:
        token = os.environ["ZENODO_TOKEN"]
    manager = ZenodoManager(token, args.sandbox)
    # Find existing depositions.
    # Zenodo's search feature is buggy, so we have to search manually.
    # repos = manager.list_depositions(
    #     q="{}:'{}'".format(*
    #         map(
    #             elasticsearch_escape,
    #             [
    #                 "metadata.custom.code:codeRepository",
    #                 os.environ["GITHUB_URL"]
    #             ]
    #         )
    #     )
    # )
    code_repo_path = parse("metadata.custom['code:codeRepository']")
    repos = [
        d for d in manager.list_depositions() if jsonpath_find_one_if_exists(
            code_repo_path,
            d
        ) == os.environ["GITHUB_URL"]
    ]
    with open(os.environ["ZENODO_JSON"]) as zenodo_json_file:
        custom_zenodo_json = json.load(zenodo_json_file)
    zenodo_data = custom_zenodo_json["zenodo"]
    title = None
    try:
        title = zenodo_data["title"]
    except KeyError:
        try:
            title = zenodo_data["metadata"]["title"]
        except KeyError:
            pass
    if title is None:
        title = os.environ["REPO_NAME"]
    zenodo_data.setdefault("metadata", {})
    zenodo_data["metadata"]["title"] = "{}: {}".format(
        title,
        os.environ["VERSION"]
    )
    #zenodo_data["metadata"]["title"] = zenodo_data["title"]
    published_at = datetime.datetime.fromisoformat(os.environ["PUBLISHED_AT"])
    zenodo_data["metadata"]["publication_date"] = published_at.strftime(
        "%Y-%m-%d"
    )
    zenodo_data["metadata"]["description"] = os.environ["RELEASE_DESCRIPTION"]
    zenodo_data["metadata"].setdefault("related_identifiers", []).append(
        {
            "identifier": os.environ["GITHUB_TAG_URL"],
            "relation": "isSupplementTo",
            "resource_type": "software",
            "scheme": "url"
        }
    )
    zenodo_data["metadata"]["version"] = os.environ["VERSION"]
    zenodo_data["metadata"].setdefault(
        "custom",
        {}
    )["code:codeRepository"] = os.environ["GITHUB_URL"]
    zenodo_data["metadata"].setdefault("imprint_publisher", "Zenodo")
    zenodo_data["metadata"]["upload_type"] = "software"
    #embed()
    if repos:
        print("Making new release.")
        new_release_id = manager.new_deposition_version(repos[0]["id"]).json()[
            "links"
        ]["latest_draft"].split("/")[-1]
        manager.update_deposition(new_release_id, zenodo_data)
        resp = manager.get_deposition(new_release_id)
    else:
        # New deposition.
        print("New deposition.")
        resp = manager.create_deposition(zenodo_data)
    if args.clean:
        for f in manager.list_files(resp.json()["id"]).json():
            manager.delete_file(resp.json()["id"], f["id"])
    bucket = resp.json()["links"]["bucket"]
    for f in args.files:
        manager.upload_file(bucket, f)
    if args.publish:
        manager.publish_deposition(resp.json()["id"])

if __name__ == "__main__":
    main()
