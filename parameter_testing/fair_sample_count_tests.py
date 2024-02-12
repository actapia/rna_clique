import argparse
import random
import itertools
import subprocess
import time
import shutil
import tqdm
from pathlib import Path
from datetime import datetime, timedelta
from dataclasses import dataclass
from queue import PriorityQueue

from find_homologs import eprint


smallest_size = 4

def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-O", "--out-dir", type=Path, default=Path("."))
    parser.add_argument("-C", "--aggressive-cleanup", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--dry-run", action="store_true")
    # parser.add_argument(
    #     "-o",
    #     "--out",
    #     type=Path,
    #     default="top_n_tests_{}_sampels"
    # )
    return parse_args_and_remainder(parser)


def parse_args_and_remainder(parser):
    args, remainder = parser.parse_known_args()
    args.fasta = []
    rem = argparse.Namespace()
    while remainder:
        if remainder[0].startswith("-"):
            key = remainder.pop(0)
            value = True
            if remainder and not remainder[0].startswith("-"):
                value = remainder.pop(0)
            setattr(rem, key, value)
        else:
            args.fasta.append(remainder.pop(0))
    return args, rem

class TopnTest():
    def __init__(self, files, out_dir, args=None):
        if args is None:
            args = []
        self.files = list(files)
        self.out_dir = out_dir
        self.args = list(args)

    def _build_command(self):
        return [
            "bash",
            "do_top_n_tests.sh"
        ] + self.args + [
            "-O",
            str(self.out_dir),            
        ] + self.files

    def run(self):
        command = self._build_command()
        proc = subprocess.Popen(
            command
        )
        start = time.perf_counter()
        proc.communicate()
        end = time.perf_counter()
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(
                proc.returncode,
                command
            )
        return end - start
            

@dataclass(order=True)
class SampleSizeElement:
    time : int = 0
    size : int = -1
    shuffles : int = 0

def make_arg_arr(d):
    args = []
    for k, v in d.__dict__.items():
        args.append(k)
        if v is not True:
            args.append(v)
    return args

class CustomRemainingBar(tqdm.tqdm):
    def __init__(
            self,
            *args,
            base_format="{l_bar}{bar}| {n_fmt}/{total_fmt}",
            **kwargs
    ):
        super().__init__(
            *args,
            bar_format=base_format + " [{elapsed}<{remaining}]",
            **kwargs
        )
        self.base_format = base_format


    def update(self, n=1):
        super().update(n)
        self.bar_format = self.base_format + " [{elapsed}<{remaining}]"
        
    def set_remaining(self, remaining):
        remaining_str = tqdm.tqdm.format_interval(remaining)
        self.bar_format = self.base_format + " [{elapsed}" + \
            f"<{remaining_str}]"

    def set_unknown(self):
        self.bar_format = self.base_format + " [{elapsed}<?]"

#from IPython import embed

def test_done(test_dir):
    return test_dir.is_dir() and all(
        out_dir.is_dir() and (out_dir / "graph.pkl").exists()
        for out_dir in test_dir.glob("out_*")
    )
        
def main():
    args, remainder = handle_arguments()
    args.out_dir.mkdir(exist_ok=True)
    if args.seed is not None:
        random.seed(args.seed)
        with open(args.out_dir / "seed", "w") as seed_file:
            seed_file.write(str(args.seed) + "\n")
    extra_args = make_arg_arr(remainder)
    subsets = {
        s: list(itertools.combinations(args.fasta, s))
        for s in range(smallest_size, len(args.fasta) + 1)
    }
    total = sum(len(s) for s in subsets.values())
    for s in subsets.values():
        random.shuffle(s)
    sample_size_queue = PriorityQueue()
    sample_sizes = []
    # from IPython import embed
    # embed()
    for i in range(smallest_size, len(args.fasta) + 1):
        elem = SampleSizeElement(0, i)
        sample_sizes.append(elem)
        sample_size_queue.put(elem)
    with CustomRemainingBar(total=total) as prog:
        while sample_size_queue:
            elem = sample_size_queue.get()
            files = subsets[elem.size].pop()
            # Do work.
            out_dir = args.out_dir / \
                f"out_{elem.size}_samples" / \
                f"shuffle_{elem.shuffles}"
            if not (args.resume and test_done(out_dir)):
                if not args.dry_run:
                    out_dir.mkdir(parents=True, exist_ok=True)
                    task = TopnTest(files, out_dir, extra_args)
                    elem.time += task.run()
                else:
                    eprint(f"Would generate {out_dir}.")
            elif args.dry_run:
                eprint(f"{out_dir} is done already. Skipping.")
            elem.shuffles += 1
            if args.aggressive_cleanup:
                if not args.dry_run:
                    shutil.rmtree(out_dir)
                else:
                    eprint(f"Would remove {out_dir}")
            # Add back to queue.
            if subsets[elem.size]:
                sample_size_queue.put(elem)
            # Compute remaining time.
            if all(el.shuffles > 0 for el in sample_sizes):
                prog.set_remaining(
                    round(
                        sum(
                            len(subsets[el.size])*el.time/el.shuffles
                            for el in sample_sizes if el.shuffles > 0
                        )
                    )
                )
            else:
                prog.set_unknown()
            prog.update()

if __name__ == "__main__":
    main()
