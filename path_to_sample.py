import re
from pathlib import Path

sample_re = re.compile(r"(.*)_top\.fasta")

def path_to_sample(x: Path | str) -> str:
    """Get a sample name from the path to a top n genes FASTA file."""
    res = sample_re.search(Path(x).name)
    if not res:
        raise PathToSampleError(
            f"Unable to parse sample name from path {x} using regular "
            f"expression {sample_re.pattern}.",
            x
        )
    return res.group(1)

def dict_path_to_sample(d):
    def inner(p):
        try:
            return d[Path(p)]
        except KeyError:
            raise PathToSampleError(
                f"Cannot get sample name for path {p}.",
                p
            )
    return inner

class PathToSampleError(Exception):
    def __init__(self, message, path):
        super().__init__(message)
        self.path = path
