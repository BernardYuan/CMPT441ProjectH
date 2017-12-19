from Bio import AlignIO
import dendropy
import numpy as np
import numpy.random as random
import scipy.optimize as opt
import copy
import sys
sys.path.append("../src")
import likelihood as llk
import tree_utils as treeutil
import sim_util as simutil
import util
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

msa = AlignIO.read("../data/fabricated_msa.fasta", "fasta")
# indel rates of each segments
segRates = [(0.4, 0.02), (80.0, 4.0)]
segRatesProbs = [0.5, 0.5]
segmentLengs = [18, 124]
msaWhole = {}
for record in msa :
    msaWhole[record.id] = copy.deepcopy(str(record.seq))

msaSep = {}
for record in msa :
    msaSep[record.id] = []
    seq = copy.deepcopy(str(record.seq))
    for l in segmentLengs :
        msaSep[record.id].append(seq[0:l])
        seq = seq[l:]
seqNames = copy.deepcopy(msaSep.keys())
alignment = copy.deepcopy(msaSep.values())

rho = 0.05
# The true tree
trueTreeStr ='[&R] (((seq1:0.1,seq2:0.1):0.1,(seq3:0.1,seq4:0.1):0.1):0.1,((seq5:0.1,seq6:0.1):0.1,(seq7:0.1,seq8:0.1):0.1):0.1);'

qMat = np.matrix([[-1.17624188,  0.54314464,  0.31904592,  0.31405132],
 [ 0.5654466,  -2.07125981,  0.62675255,  0.87906067],
 [ 0.15170778,  0.28626922, -1.05652458,  0.61854757],
 [ 0.08300798,  0.22318344,  0.34382518, -0.65001661]])

distMat = []
for i in range(len(seqNames)) :
    distMat.append([])
    for j in range(i) :
        res = opt.minimize_scalar(llk.GeoPIP_likelihood_2list, args = (alignment[i], alignment[j], segRates, qMat))
        # print seqNames[i], "->", seqNames[j], ":", res.x
        distMat[i].append(res.x)
    distMat[i].append(0)

distanceMatrix = DistanceMatrix(seqNames, distMat)
print distanceMatrix

# DistanceTreeConstructor().nj()
trueTree = dendropy.Tree.get(data = trueTreeStr, schema="newick")
truenllk = llk.GeoPIP_likelihood(msaWhole.values(), msaWhole.keys(), segmentLengs, segRates, segRatesProbs, trueTree, qMat, rho)
print "Tree Tree:", truenllk

njTreeStr = treeutil.build_tree_UPGMA(msa, distanceMatrix)
print njTreeStr

treeSet = treeutil.NNITreeSurgery(njTreeStr)
treeSet.insert(0, njTreeStr)
for ts in treeSet :
    tT = dendropy.Tree.get(data=ts, schema="newick")
    nllk = llk.GeoPIP_likelihood(msaWhole.values(), msaWhole.keys(), segmentLengs, segRates, segRatesProbs, tT, qMat, rho)
    print nllk
