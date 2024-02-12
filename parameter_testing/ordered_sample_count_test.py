import argparse
from pathlib import Path
from fair_sample_count_tests import (
    parse_args_and_remainder,
    TopnTest,
    smallest_size
)

def handle_arguments(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser()
    parser.add_argument("-O", "--out-dir", type=Path, default=Path("."))
    return parse_args_and_remainder(parser)

def prefixes(l, min_size=0):
    for i in range(min_size, len(l) + 1):
        yield l[:i]

class OrderedSampleCountTest:
    def __init__(self, files, out_dir, min_size=smallest_size, args=None):
        if args is None:
            args = []
        self.files = list(files)
        self.out_dir = out_dir
        self.args = list(args)
        self.min_size = min_size

    def run(self):
        for files in prefixes(self.files, self.min_size):
            out_dir = self.out_dir / "out_{}_samples".format(len(files))
            task = TopnTest(files, out_dir, self.args)
            task.run()
