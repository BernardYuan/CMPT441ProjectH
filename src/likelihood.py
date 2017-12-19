# external modules
import numpy as np
import scipy as sp
import scipy.special as special
import scipy.stats as stats
import copy
# Local modules
import util
#testing variables
# msa = [['A','T','G','C'],['A','-','A','-'],['-','-','C','C']]
# seglen = [1,1,2]
chIndex = {'A':0, 'C':1, 'G':2, 'T':3, '-':4}
def pc_from_tree_with_pv_and_fv(treeWithPvFv) :
	pc = 0
	for node in treeWithPvFv.nodes() :
		pc = pc + node.pv * node.fv
	return pc

# Use Felsenstein peel algorithm to compute ~f
# vector and ~f at each inner node
def felsenstein_peel(qMatExt, piExt, ftv1, bl1, ftv2, bl2) :
	pMat1 = sp.linalg.expm(bl1 * qMatExt)
	pMat2 = sp.linalg.expm(bl2 * qMatExt)
	ftv = np.dot(pMat1, ftv1) * np.dot(pMat2, ftv2)
	ft = np.dot(piExt, ftv)
	return ftv, ft

# Compute ~f at leaf nodes
def ftilde_leaf(piExt, ch) :
	idx = chIndex[ch]
	ftv = np.zeros(5)
	ftv[idx] = 1
	ft = piExt[idx]
	return ftv, ft

# Compute the term fv by col
def fv_by_site(msaSite, seqNames, treeWithBeta, qMatExt, piExt, dRate) :
	nSeq = len(msaSite)
	isCphi = (msaSite.count('-') == len(msaSite))
	# make the character dictionary
	chDict = dict(zip(seqNames, msaSite))

	# If not Cphi, find the common ancestors and
	# the path from root to the most recent common ancestors
	if not isCphi:
		# Set S in the support information
		nonEmptySeqNames = [seqNames[i] for i in range(nSeq) if msaSite[i]!='-']
		# The youngest common ancestor of S
		youngestCommonAncestor = treeWithBeta.mrca(taxon_labels=nonEmptySeqNames)
		# Set A in the support information
		nonZeroFvNodes = util.node_path(treeWithBeta, youngestCommonAncestor)
	
	for node in treeWithBeta.postorder_node_iter() :
		if node.is_leaf() :
			nodeName = node._get_node_token()
			ch = chDict[nodeName]

			# ~f vector and sum over all neucleotides
			ftilde_v, ftilde = ftilde_leaf(piExt, ch)
			if isCphi :
				fv = 1 + node.beta * (ftilde - 1)
			elif node in nonZeroFvNodes :
				fv = ftilde * node.beta
			else :
				fv = 0
			node.ftilde_v = ftilde_v
			node.ftilde = ftilde
			node.fv = fv
		else :
			children = node.child_nodes()
			child1 = children[0]
			child2 = children[1]

			ftv1 = child1.ftilde_v
			ftv2 = child2.ftilde_v
			bl1 = child1.edge_length
			bl2 = child2.edge_length

			ftv, ft = felsenstein_peel(qMatExt, piExt, ftv1, bl1, ftv2, bl2)
			# Compute fv at internal points
			if isCphi :
				fv = 1 + node.beta * (ft - 1)
			elif node in nonZeroFvNodes :
				fv = node.beta * ft
			else :
				fv = 0
			node.ftilde_v = ftv
			node.ftilde = ft
			node.fv = fv

# Compute the prior of insertion and survival rate in the tree
def pv_and_beta_recursion_tree(tree, dRate) :
	tau = tree.length()
	dRateInverse = 1.0/dRate
	nuBar = 1.0 / (tau + dRateInverse)
	#interator of the tree
	iterTree = tree.preorder_node_iter()
	# Work on the root of the tree
	root = iterTree.next()
	index = 0
	root.index = index
	root.pv = nuBar * dRateInverse
	root.beta = 1

	for node in iterTree:
		index = index + 1
		node.index = index

		node.pv = node.edge_length * nuBar
		betaMu = node.edge_length * dRate
		if betaMu <= 0:
			node.beta = 0
		else :
			node.beta = 1. / betaMu * (1 - np.exp(-betaMu))

# Compute the likelihood of each site
def PIP_likelihood_single_site(msaSite, seqNames, tree, qMatExt, piExt, dRate) :
	fv_by_site(msaSite, seqNames, tree, qMatExt, piExt, dRate) 
	pc = pc_from_tree_with_pv_and_fv(tree)
	return pc

# Compute the likelihood of multiple sites
def PIP_likelihood_multiple_sites(msaBySites, seqNames, indelRate, tree, qMatrix) :
	iRate = indelRate[0]
	dRate = indelRate[1]
	piProb = util.pi_from_qmat(qMatrix)
	piProbExt = np.append(piProb, 0)
	qMatExt = util.q_to_qext(qMatrix, dRate)

	# Compute pv and beta in the tree
	pv_and_beta_recursion_tree(tree, dRate)

	# Compute PIP(c) on one segment
	pcMultiSite = 0
	for site in msaBySites:
		pcOneSite = PIP_likelihood_single_site(site, seqNames, tree, qMatExt, piProbExt, dRate)
		pcMultiSite = pcMultiSite + np.log(pcOneSite)
	
	# Compute psi 
	tau = tree.length()
	nu = iRate * (tau + 1.0/dRate)
	nSeqs = len(seqNames)
	nCols = len(msaBySites)
	cphi = []
	for i in range(nSeqs) :
		cphi.append('-')
	# Compute P(C_0)
	pc0 = PIP_likelihood_single_site(cphi, seqNames, tree, qMatExt, piProbExt, dRate)

	# Compute log Psi
	logPsi = -np.sum(np.log(np.arange(1, nCols + 1))) + np.log(nu) * nCols + (pc0 - 1) * nu
	return logPsi + pcMultiSite

# Compute P_{tau}(m) over all segments
def GeoPIP_likelihood(msa, seqNames, segLength, segRates, segRateProbs, tree, qMatrix, rho) :
	msaBySeg = util.divide_MSA_by_segment(msa, segLength)

	#Initial likelihood
	llh = 0

	#Compute Likelihood by segment
	# Multiply them together with omega
	segIDs = sorted(msaBySeg.keys())
	for segid in segIDs :
		segLLH = PIP_likelihood_multiple_sites(msaBySeg[segid], seqNames, segRates[segid], tree, qMatrix) 
		llh = llh - segLLH
	
	# add log of omegas
	llh = llh - np.sum(np.log(segRateProbs))

	# |Beta|
	nSegs = len(segIDs)
	# take into account the term with geometric parameter
	logGeo = (nSegs-1)*np.log(rho) + np.log(1 - rho)
	llh = llh - logGeo
	return llh

def pv_prior(b1, b2, dRate) :
	tau = b1 + b2
	denom = tau + 1./dRate
	p1 = b1/denom
	p2 = b2/denom
	p0 = (1./dRate) / denom
	return p1, p2, p0

def fv(c1, c2, piExt, pMatExt, b, dRate, cListExt) :
	index1 = cListExt.index(c1)
	index2 = cListExt.index(c2)
	pic1 = piExt[index1]
	pic2 = piExt[index2]
	logic1 = (c1 != '-')
	logic2 = (c2 != '-')
	if logic1:
		ft0 = pic1 * pMatExt[index1, index2]
	else:
		ft0 = 0
	f0 = ft0
	beta = (1 - np.exp(-dRate * b)) / (b * dRate)
	if ((not logic1) and logic2):
		f2 = beta * pic2
	elif ((not logic1) and (not logic2)):
		f2 = 1 - beta
	else:
		f2 = 0
	if (logic1 or logic2):   # not really necesarry
		f1 = 0   # not really necesarry
	else:   # not really necesarry
		f1 = 1   # not really necesarry
	return f1, f2, f0

def prob_c_all(piExt, pMat, b, dRate, cListExt) :
	pc = {}
	pvPrior = pv_prior(0, b, dRate)
	for c1 in cListExt :
		pc[c1] = {}
		for c2 in cListExt :
			fvProb = fv(c1, c2, piExt, pMat, b, dRate, cListExt)
			pc[c1][c2] = np.dot(pvPrior, fvProb)
	return pc

# Compute likelihood between two sequences
def logprob_m(seq1, seq2, b, pcAll, iRate, dRate) :
	tau = b
	nu = iRate * (tau + 1./dRate)
	nSites = len(seq1)
	if nSites == 0 :
		return (pcAll['-']['-'] - 1) * nu
	else :
		logphi = -np.sum(np.log(np.arange(1, nSites+1))) + nSites * np.log(nu) + (pcAll['-']['-']-1) * nu
		logpc = 0
		for i in range(nSites) :
			logpc += np.log(pcAll[seq1[i]][seq2[i]])
		return logphi + logpc

def PIP_likelihood_2list(b, seq1, seq2, indelRate, qmat, cList=['A', 'C', 'G', 'T']) :
	# probability of insertion
	insRate = indelRate[0]
	delRate = indelRate[1]
	qext = util.q_to_qext(qmat, delRate)
	pmat = sp.linalg.expm(b * qext)
	pi = util.pi_from_qmat(qmat)
	piExt = np.append(pi, [0])
	cListExt = cList + ['-']
	pcAll = prob_c_all(piExt, pmat, b, delRate, cListExt)
	return logprob_m(seq1, seq2, b, pcAll, insRate, delRate)

# Compute the likelihood between two sequences with GeoPIP model
# b: integer
# seq1: [[seg1], [seg2], ...]
# seq2: [[seg1], [seg2], ...]
# segRates: [(ins1, del1), (ins2, del2)]
# qmat: the rate matrix
def GeoPIP_likelihood_2list(b, seq1, seq2, segRates, qMat) :
	if b <= 0:    # set lower bound, do not accept negative values
		return 1.e10
	if b > 100:    # set upper bound, do accept bigger values, but truncate
		return 1.e-10    # if bigger than 100, then return 100
	nllk = 0
	nSegs = len(segRates)
	for i in range(nSegs) :
		seg1 = copy.deepcopy(seq1[i])
		seg2 = copy.deepcopy(seq2[i])
		indelRate = segRates[i]
		nllk += -PIP_likelihood_2list(b, seg1, seg2, indelRate, qMat)
	return nllk