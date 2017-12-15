import numpy as np

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