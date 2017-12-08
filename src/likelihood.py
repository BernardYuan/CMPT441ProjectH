# external modules
import numpy as np
import scipy as sp
import scipy.special as special

# Local modules
import util

#testing variables
# msa = [['A','T','G','C'],['A','-','A','-'],['-','-','C','C']]
# seglen = [1,1,2]

chIndex = {'A':0, 'T':1, 'G':2, 'C':3, '-':4}

def pc_from_tree_with_pv_and_fv(treeWithPvFv) :
	pc = 0
	for node in treeWithPvFv :
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
	nSeq = len(msa)
	isCphi = (msa.count('-') == len(msaSite))
	# make the character dictionary
	chDict = dict(zip(seqNames, msaSite))

	# If not Cphi, find the common ancestors and
	# the path from root to the most recent common ancestors
	if not isCphi:
		# Set S in the support information
		nonEmptySeqNames = [segNames[i] for i in range(nSeq) if msaSite[i]!='-']
		# The youngest common ancestor of S
		youngestCommonAncestor = treeWithBeta.mrca(taxon_labels=nonEmptySeqNames)
		# Set A in the support information
		nonZeroFvNodes = util.nodePath(treeWithBeta, youngestCommonAncestor)
	
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
			child1, child2 = node.child_nodes()
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

# Compute psi provided z and k
def compute_Psi(z, k, nu) :
	# first term in psi
	k_facorial = special.factorial(k)
	inverse_k_factorial = 1.0 / k_factorial
	# second term in psi
	nu_power_k = np.power(nu, k)
	# third term in psi
	exp_zminus1_nu = np.exp((z-1)*nu)

	return inverse_k_factorial * nu_power_k * exp_zminus1_nu


# Compute the likelihood of multiple sites
def PIP_likelihood_multiple_sites(msaBySites, seqNames, indelRate, tree, qMatrix) :
	iRate = indelRate[0]
	dRate = indelRate[1]
	piProb = util.pi_from_qmat(qMatrix)
	piProbExt = np.append(piProb, 0)
	qMatExt = util.q_to_qext(qMatrix, dRate)

	# Compute pv and beta in the tree
	pv_and_beta_recursion_tree(tree, dRate)

	# Compute PIp(c) on one segment
	pcMultiSite = 1
	for site in msaBySites:
		pcMultiSite = pcMultiSite * PIP_likelihood_single_site(site, seqNames, tree, qMatExt, piProbExt, dRate)
	
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
	psi = compute_Psi(pc0, nCols, nu)

	return psi * pcMultiSite

# Compute P_{tau}(m) over all segments
def GeoPIP_likelihood(msa, seqNames, segLength, segRates, segRateProbs, tree, qMatrix, rho) :
	msaBySeg = util.divide_MSA_by_segment(msa, segLength)
	#Initial likelihood
	llh = 1
	#Compute Likelihood by segment
	# Multiply them together with omega
	segIDs = sorted(msaBySeg.keys())
	for segid in segIDs :
		segLLH = PIP_likelihood_by_Sites(msaBySeg[segid], segRates[segid], tree, qMatrix) 
		llh = llh * segRateProbs[segid] * segLLH

	# |Beta|
	nSegs = len(segIDs)

	# multiply the term with geometric parameter
	llh = llh * np.power(1-rho, nSegs-1) * rho
	return llh

