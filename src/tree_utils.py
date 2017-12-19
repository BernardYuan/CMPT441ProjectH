from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import NNITreeSearcher
from Bio.Phylo.TreeConstruction import Scorer
from Bio import Phylo
from cStringIO import StringIO
import util
import numpy as np
import numpy.random as random
from dendropy import Tree

JukesCantorQmatrix = [
    [-3.0, 1.0, 1.0, 1.0],
    [1.0, -3.0, 1.0, 1.0],
    [1.0, 1.0, -3.0, 1.0],
    [1.0, 1.0, 1.0, -3.0]]

treeTest = "[&R] (((seq3:0.15533,seq4:0.03316)Inner4:0.08312,(seq2:0.09434,seq1:0.13294)Inner5:0.06637)Inner6:0.01930,((seq5:0.05086,seq6:0.08883)Inner2:0.05831,(seq7:0.07312,seq8:0.07847)Inner1:0.10408)Inner3:0.07520):0.00000;"

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
                    if row.seq[i]!='-' and col.seq[i]!= '-' and row.seq[i] != col.seq[i] :
                        diff = diff + 1
                JDdist = -0.75 * np.log(1 - 4./3 * (1.0 * diff / strLen))
                matrix[rowIdx].append(JDdist)
        rowIdx = rowIdx + 1
    return DistanceMatrix(names, matrix)
def build_tree_UPGMA(msa, distanceMatrix=None) :
    if not distanceMatrix :
        distCalculator = DistanceCalculator("identity")
        distanceMatrix = distCalculator.get_distance(msa)
    # Construct the tree with the distance Matrix
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distanceMatrix)
    # Make the tree rooted
    #tree.root_at_midpoint()
    #return newick format
    return "[&R] " + tree.format("newick").strip()

def build_tree_NJ(msa, distanceMatrix=None) :
    if not distanceMatrix :
        distCalculator = DistanceCalculator("identity")
        distanceMatrix = distCalculator.get_distance(msa)
    # Construct the tree with the distance Matrix
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distanceMatrix)
    # Make the tree rooted
    #tree.root_at_midpoint()
    #return newick format
    return "[&R] " + tree.format("newick").strip()

def NNITreeSurgery(treeStr, schema="newick") :
    tree = Phylo.read(StringIO(treeStr), schema)

    nniTreeSearcher = NNITreeSearcher(Scorer())
    treeList = nniTreeSearcher._get_neighbors(tree)
    allTreeList = []
    for t in treeList :
        tstr = t.format("newick")
        tstr.strip()
        allTreeList.append("[&R] " + tstr[0:len(tstr)-6] + ";")
    return allTreeList

def node_distance(offspring, ancestor) :
    node = offspring
    dist = 0
    while node != ancestor :
        dist += node.edge_length
        node = node.parent_node
    return dist

def leaf_distance_from_tree(tree) :
    leavesName = []
    for leaf in tree.leaf_nodes() :
        leavesName.append(leaf.taxon.label)
    distMatrix = []
    nLeaf = len(leavesName)
    for i in range(nLeaf) :
        distMatrix.append([])
        for j in range(i) :
            name1 = leavesName[i]
            name2 = leavesName[j]
            node1 = tree.find_node_with_taxon_label(name1)
            node2 = tree.find_node_with_taxon_label(name2)
            ancestor = tree.mrca(taxon_labels=[name1, name2])
            dist1 = node_distance(node1, ancestor)
            dist2 = node_distance(node2, ancestor)
            distMatrix[i].append(dist1 + dist2)
        distMatrix[i].append(0)
    return DistanceMatrix(leavesName, distMatrix)

def generate_balanced_tree_with_fixed_edge_length(numLeaf, bLen=0.5, prefix="seq") :
    leafNames = [prefix+str(i) for i in range(numLeaf)]
    nodeNames = [name for name in leafNames]

    while len(nodeNames) > 1 :
        node1 = nodeNames[0]
        node2 = nodeNames[1]

        nodeNew = '(' + node1 + ':' + str(bLen) + ',' + node2 + ':' + str(bLen) + ')'
        nodeNames = nodeNames[2:]
        nodeNames.append(nodeNew)
    
    return '[&R] ' + nodeNames[0] + ';'