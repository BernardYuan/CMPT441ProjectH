from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio import Phylo
import numpy as np

def JukesCantorDistanceMatrix(msa) :
    names = [seq.id for seq in msa]
    matrix = []

    rowIdx = 0
    for row in msa:
        matrix.append([])
        for col in msa:
            if col.id == row.id :
                matrix[rowIdx].append(0)
                break
            else :
                strLen = len(row.seq)
                diff = 0
                for i in range(strLen) :
                    if row.seq[i] != col.seq[i] :
                        diff = diff + 1
                JDdist = -0.75 * np.log(1 - 4./3 * (1.0 * diff / strLen))
                matrix[rowIdx].append(JDdist)
        rowIdx = rowIdx + 1
    return DistanceMatrix(names, matrix)

def build_tree_NJ(msa, distanceMatrix=None) :
    if not distanceMatrix :
        distCalculator = DistanceCalculator("identity")
        distanceMatrix = distCalculator.get_distance(msa)

    # Construct the tree with the distance Matrix
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distanceMatrix)

    # Make the tree rooted
    tree.root_at_midpoint()
    #return newick format
    return "[&R] " + tree.format("newick").strip()