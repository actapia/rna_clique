import re
from pathlib import Path

sample_re = re.compile(r"(.*)_top\.fasta")

def path_to_sample(x):
    return sample_re.search(Path(x).name).group(1)