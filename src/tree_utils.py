from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo


def build_tree_NJ(msaFile, msaFileType, distanceModel="identity") :
    # read the msa and compute pairwise distance
    msa = AlignIO.read(msaFile, msaFileType)
    distCalculator = DistanceCalculator(distanceModel)

    # Construct the tree
    constructor = DistanceTreeConstructor(distCalculator, 'nj')
    tree = constructor.build_tree(msa)

    #return newick format
    return tree.format("newick")