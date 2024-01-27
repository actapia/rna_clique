import argparse
import re
import Bio.Phylo
import dendropy
import dendropy.calculate

from io import StringIO
from pathlib import Path
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from dendropy.calculate.treecompare import symmetric_difference
from tqdm import tqdm
from filtered_distance import SampleSimilarity
from phylo_utils import tril_jagged

def phylo_to_dendropy(tree, ns=None):
    sio = StringIO()
    Bio.Phylo.write(tree, sio, format="newick")
    return dendropy.Tree.get(data=sio.getvalue(), schema="newick", taxon_namespace=ns)

def multi_phylo_to_dendropy(*trees, ns=None):
    ns = None
    for tree in trees:
        new = phylo_to_dendropy(tree, ns=ns)
        ns = new.taxon_namespace
        yield new


def handle_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("analysis_root", type=Path)
    parser.add_argument("ground_truth", type=Path)
    return parser.parse_args()

def main():
    args = handle_arguments()
    graph_path = args.analysis_root / "graph.pkl"
    comparisons = list((args.analysis_root / "od2").glob("*.pkl"))
    sim = SampleSimilarity.from_filenames(graph_path, tqdm(comparisons))
    dis_matrix = sim.get_dissimilarity_matrix()
    sim_name_re = re.compile("(T.*)_top")
    sample_taxa = [
        sim_name_re.search(Path(s).stem).group(1) for s in sim.samples
    ]
    phylo_dis_mat = DistanceMatrix(sample_taxa, tril_jagged(dis_matrix))
    constructor = DistanceTreeConstructor()
    ground_truth = next(Bio.Phylo.parse(args.ground_truth, format="newick"))
    nj_tree = constructor.nj(phylo_dis_mat)
    assert symmetric_difference(
        *multi_phylo_to_dendropy(nj_tree, ground_truth)
    ) == 0

if __name__ == "__main__":
    main()
