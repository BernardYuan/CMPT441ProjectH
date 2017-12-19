from Bio import AlignIO
import dendropy
import numpy as np
import numpy.random as random
import sys

sys.path.append("../src")
import likelihood as llk
import tree_utils as treeutil
import util

# read in the synthetic data
msa = AlignIO.read("../data/fabricated_msa.fasta", "fasta")
# indel rates of each segments
segRates = [(0.4, 0.02), (80.0, 4.0)]
segRatesProbs = [0.5, 0.5]
segmentLengs = [18, 124]
seqNames = [record.id for record in msa]
alignment = [record.seq for record in msa]
rho = 0.05
# The true tree
treeStr ='[&R] (((seq1:0.1,seq2:0.1):0.1,(seq3:0.1,seq4:0.1):0.1):0.1,((seq5:0.1,seq6:0.1):0.1,(seq7:0.1,seq8:0.1):0.1):0.1);'
qMat = np.matrix([[-1.17624188,  0.54314464,  0.31904592,  0.31405132],
 [ 0.5654466,  -2.07125981,  0.62675255,  0.87906067],
 [ 0.15170778,  0.28626922, -1.05652458,  0.61854757],
 [ 0.08300798,  0.22318344,  0.34382518, -0.65001661]])

# Build the distance matrix and use the distance
# matrix to construct a tree with NJ
# distanceMatrix = treeutil.JukesCantorDistanceMatrix(msa)
# treeStr = treeutil.build_tree_NJ(msa, distanceMatrix)

# All possible indelRates
#indelRates = [(80.0, 2.0), (4.0, 0.1)]
#indelRateProbs = [0.5, 0.5]

# use a Jukes and Cantor matrix for qMat
#qMat = treeutil.JukesCantorQmatrix
# number of segments
#m = 5
# segmentLengs = util.generate_segment_lengths(m, len(alignment[0]))
# segRates, segRatesProbs = util.generate_segment_indel_rates(m, indelRates, indelRateProbs)

# from the tree string and generate a set of trees with NNI
treeSetNNI = treeutil.NNITreeSurgery(treeStr)
# Base tree nllk
tree = dendropy.Tree.get(data=treeStr, schema="newick")
nllk = llk.GeoPIP_likelihood(alignment, seqNames, segmentLengs, segRates, segRatesProbs, tree, qMat, rho)
print "True tree:", nllk, "  tree string:", treeStr
# surgery tree
index = 0
for treestr in treeSetNNI :
    index = index + 1
    tree = dendropy.Tree.get(data=treestr, schema="newick")
    nllk = llk.GeoPIP_likelihood(alignment, seqNames, segmentLengs, segRates, segRatesProbs, tree, qMat, rho)
    print "tree", index, ":", nllk, "  tree string:", treestr