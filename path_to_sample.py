import re
from pathlib import Path

sample_re = re.compile(r"(.*)_top\.fasta")

def path_to_sample(x: Path | str) -> str:
    """Get a sample name from the path to a top n genes FASTA file."""
    return sample_re.search(Path(x).name).group(1)
