import numpy as np
import numpy.random as random

# Divide a MSA by segment
# Input: 
#	msa: a multiple sequence alignment
#	segLength: the length of each segment
# Output:
#	a map where key is segment ID and value is a MSA segment.
#Example:
#Input:
# 	msa: 
#		[['ATGC'],
#		 ['A-A-'],
#		 ['--CC']]
#	segLength:
#		[1,2,1]
#Output:
#	{
#	0:[['A', 'A', '-']],
#	1:[['T', '-', '-'],
#	   ['G', 'A', 'C']],
#   2:[['C', '-', 'C']]
#	}
def divide_MSA_by_segment(msa, segLength) :
	msa_by_seg = {}
	segID = 0
	colBase = 0

	for seglen in segLength :
		msaSeg = []
		for i in range(seglen) :
			colIdx = colBase + i;
			segCol = []
			for seq in msa :
				segCol.append(seq[colIdx]);
			msaSeg.append(segCol)

		colBase = colBase + seglen
		msa_by_seg[segID] = msaSeg
		segID = segID + 1
	
	return msa_by_seg

# q to q_ext, add deletion rate to it
def q_to_qext(qMat, dRate):
	m,n = np.shape(qMat)
	qMatExt = np.zeros((m+1, n+1), dtype=float)
	qMatExt[:n, :n] = qMat
	qMatExt[:n, n] = dRate

	np.fill_diagonal(qMatExt, 0)
	rowSum = qMatExt.sum(axis=1)
	np.fill_diagonal(qMatExt, -rowSum)

	return qMatExt

# Compute the stationary distribution pi from qMat
def pi_from_qmat(qMat) :
	m,n = np.shape(qMat)

	qMatExt = np.hstack([qMat, np.ones((m, 1))])
	b = np.zeros((m+1, ))
	b[-1] = 1

	res = np.linalg.lstsq(qMatExt.T, b)
	x = res[0]
	return x

# root to node path
def node_path(tree, node) :
	nodePath = []
	nodePtr = node
	while node is not tree.seed_node :
		nodePath.append(node)
		node = node.parent_node
	nodePath.append(node)

	return nodePath

# Generate segment lengths
def generate_segment_lengths(m, sequenceLength) :
	# Randomly generate segments
	segSites = random.uniform(0,1,m-1)
	segSites.sort()
	segSites = np.append(segSites, [1])
	segSites = np.floor(sequenceLength * segSites)

	segmentLengs = []
	for index in range(1, len(segSites)) :
		curLength = segSites[index] - segSites[index-1]
		segmentLengs.append(curLength)
	segmentLengs = [int(l) for l in segmentLengs]
	return segmentLengs

def generate_segment_indel_rates(m, indelRates, indelRateProbs) :
	indelRateProbCumSum = []
	for index in range(len(indelRateProbs)) :
		if index == 0 :
			indelRateProbCumSum.append(indelRateProbs[index])
		else :
			temSum = indelRateProbs[index] + indelRateProbCumSum[index-1]
			indelRateProbCumSum.append(temSum)
	### Allocate indel rates
	indelSeeds = random.uniform(0,1,5)
	segRates = []
	segRatesProbs = []
	for seed in indelSeeds :
		for index in range(len(indelRateProbCumSum)) :
			if seed <= indelRateProbCumSum[index] :
				segRates.append(indelRates[index])
				segRatesProbs.append(indelRateProbs[index])

	return segRates,segRatesProbs