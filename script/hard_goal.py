from Bio import AlignIO
import dendropy
import numpy as np
import numpy.random as random
import sys

sys.path.append("../src")
import likelihood as llk
import tree_utils as treeutil
import util

# read in the molluscan data
msa = AlignIO.read("../data/fabricated_msa.fasta", "fasta")
seqNames = [record.id for record in msa]
alignment = [record.seq for record in msa]
# Build the distance matrix and use the distance
# matrix to construct a tree with NJ
distanceMatrix = treeutil.JukesCantorDistanceMatrix(msa)
treeStr = treeutil.build_tree_NJ(msa, distanceMatrix)
# All possible indelRates
indelRates = [(80.0, 2.0), (4.0, 0.1)]
indelRateProbs = [0.5, 0.5]
# use a Jukes and Cantor matrix for qMat
qMat = treeutil.JukesCantorQmatrix
# number of segments
m = 5
segmentLengs = util.generate_segment_lengths(m, len(alignment[0]))
segRates, segRatesProbs = util.generate_segment_indel_rates(m, indelRates, indelRateProbs)
# Geometric parameter
rho = 0.5

# from the tree string and generate a set of trees with NNI
treeSetNNI = treeutil.NNITreeSurgery(treeStr)

# Base tree nllk
tree = dendropy.Tree.get(data=treeStr, schema="newick")
nllk = llk.GeoPIP_likelihood(alignment, seqNames, segmentLengs, segRates, segRatesProbs, tree, qMat, rho)
print "Base tree:", nllk

# surgery tree
index = 0
for treestr in treeSetNNI :
    index = index + 1
    tree = dendropy.Tree.get(data=treestr, schema="newick")
    nllk = llk.GeoPIP_likelihood(alignment, seqNames, segmentLengs, segRates, segRatesProbs, tree, qMat, rho)
    print "tree", index, ":", nllk