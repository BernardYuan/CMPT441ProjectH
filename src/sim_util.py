import numpy as np
import numpy.random as random
import time
import copy
import scipy as scp
import scipy.special as special
import bisect
import util

# generate the root node
def initial_node(m, segRates, qMat, cList = ['A', 'C', 'G', 'T']) :
    segLen = []
    for idRate in segRates :
        iRate = idRate[0]
        dRate = idRate[1]
        segLen.append(random.poisson(iRate/dRate))

    pi = util.pi_from_qmat(qMat)
    piSum = np.array(pi).cumsum()
    res = []
    for seglen in segLen :
        segprobs = random.uniform(size = seglen)
        sequence = []
        for p in segprobs :
            sequence.append((time.time(), cList[bisect.bisect(piSum, p)]))
        res.append(sequence)
    return res

def seg_evolve_withPIP(seq, blen, indelRate, qmat, cList = ['A', 'C', 'G', 'T']) :
    insertRate = indelRate[0]
    deleteRate = indelRate[1]
    subRateVec = np.array([-qmat[cList.index(pair[1]), cList.index(pair[1])] for pair in seq if pair[1] != '-'])
    pi = util.pi_from_qmat(qmat)
    piCum = np.array(pi).cumsum()

    # deep copy the original segment
    newSeq = copy.deepcopy(seq)
    # non empty column index
    seqIndexNonEmpty = [i for i in range(len(newSeq)) if newSeq[i][1]!='-']
    nonEmptyLen = len(seqIndexNonEmpty)
    # first generate insertion time
    tIns = random.exponential(scale = 1.0/insertRate)
    # if the sequence is empty
    if nonEmptyLen == 0:
        # Nothing happens if the branch length is not long enough
        # Otherwise, insertion happens 
        if blen < tIns :
            return newSeq, blen - tIns
        else :
            p = random.uniform()
            newSeq.insert(0, (time.time(), cList[bisect.bisect(piCum, p)]))
            return newSeq, blen - tIns
    # deletion and substitution rate
    tDelVec = random.exponential(scale = 1.0/deleteRate, size=nonEmptyLen)
    tSubVec = random.exponential(scale = 1.0/subRateVec)
    delMinIndex = tDelVec.argmin()
    delMin = tDelVec.min()
    subMinIndex = tSubVec.argmin()
    subMin = tSubVec.min()
    # minimum time 
    tMinVec = np.array([tIns, delMin, subMin])
    tMin = tMinVec.min()
    tMinIndex = tMinVec.argmin()
    # if time is not enough for any change
    if blen < tMin :
        return newSeq, blen - tMin
    # Insertion
    if tMinIndex == 0 :
        seqIndexNonEmpty.insert(0, seqIndexNonEmpty[0]-1)
        insPos = random.choice(seqIndexNonEmpty)
        insProb = random.uniform()
        chNew = cList[bisect.bisect(piCum, insProb)]
        newSeq.insert(insPos+1, (time.time(), chNew))
    # Deletion
    elif tMinIndex == 1 :
        delPos = seqIndexNonEmpty[delMinIndex]
        newSeq[delPos]= (newSeq[delPos][0], '-')
    # substitution
    else :
        subPos = seqIndexNonEmpty[subMinIndex]
        chToSub = newSeq[subPos][1]
        # Transitional probability
        tProb = scp.linalg.expm(tMin * qmat)[cList.index(chToSub), ]
        tProbCum = tProb.cumsum()
        subProb = random.uniform()
        chNew = cList[bisect.bisect(tProbCum, subProb)]
        newSeq[subPos] = (newSeq[subPos][0],chNew)
    return newSeq, blen - tMin

def sequence_evolve(seqAllSegs, blen, segRates, qmat, cList=['A', 'C', 'G', 'T']) :
    nSegs = len(seqAllSegs)
    for i in range(nSegs) :
        t = blen
        newSeg = copy.deepcopy(seqAllSegs[i])
        indelRate = segRates[i]
        while t > 0 :
            newSeg, t = seg_evolve_withPIP(newSeg, t, indelRate, qmat, cList=['A', 'C', 'G', 'T'])
        seqAllSegs[i] = newSeg
    return seqAllSegs

def simulate_evolution(tree, m, segRates, qMat) :
    for node in tree.preorder_node_iter() :
        # root node
        if node.parent_node is None :
            node.value = initial_node(m, segRates, qMat)
        else :
            pNode = node.parent_node
            pSeq = copy.deepcopy(pNode.value)
            blen = node.edge_length
            node.value = sequence_evolve(pSeq, blen, segRates, qMat)

def collect_leaf_alignment(tree) :
    nSegs = len(tree.seed_node.value)
    seqNames = []
    for leaf in tree.leaf_nodes() :
        seqNames.append(leaf.taxon.label)
    # Get the timestamp index of each col of each seg
    # each seg in segs is a map: seqName -> timestamp
    segs = []
    for i in range(nSegs) :
        timestamps = {}
        for leaf in tree.leaf_nodes() :
            timestamps[leaf.taxon.label] = list(zip(*(leaf.value[i]))[0])
        segs.append(timestamps)
    # Build a dictionary of each sequence,
    # which is a map: seqname -> timestamp -> neuleotide
    sequences = {}
    for leaf in tree.leaf_nodes() :
        leafValue = []
        for pairs in leaf.value :
            leafValue += copy.deepcopy(pairs)
        sequences[leaf.taxon.label] = dict(leafValue)
    # deal with each seg
    alignedTimeStamps = []
    for i in range(nSegs) :
        timeStamps = segs[i]
        alignedSegTS = []
        while True :
            # find the topologically first timestamp as the index to process
            # First find a candidate set of timestamps that can be topologically first one
            candidateTs = set()
            for name in seqNames :
                if len(timeStamps[name]) > 0 : 
                    candidateTs.add(timeStamps[name][0])
            if len(candidateTs) == 0 :
                break
            # verify each TS and see whether it is topologically first
            verifiedTs = set()
            for cT in candidateTs : 
                indexSum = 0
                for name in seqNames :
                    if len(timeStamps[name]) > 0 and (cT in timeStamps[name]):
                        indexSum += timeStamps[name].index(cT)
                if indexSum == 0:
                    verifiedTs.add(cT)
            if len(verifiedTs) == 0 :
                print ("==================Error In Verified Set===========================")
                exit()
            # add a col into the alignment of the seg
            for vT in verifiedTs :
                alignedSegTS.append(vT)
                for name in seqNames : 
                    if len(timeStamps[name]) > 0 and timeStamps[name][0] == vT :
                        timeStamps[name] = timeStamps[name][1:]
        alignedTimeStamps.append(alignedSegTS)
    
    multipleAlignment = {}
    for name in seqNames :
        multipleAlignment[name] = []
        for i in range(nSegs) :
            multipleAlignment[name].append('')
    for i in range(nSegs) :
        segTs = alignedTimeStamps[i]
        for ts in segTs :
            for name in seqNames :
                if ts in sequences[name] :
                    multipleAlignment[name][i] += sequences[name][ts]
                else :
                    multipleAlignment[name][i] += '-'
    return multipleAlignment

# remove gaps in MSA
def msa_remove_gaps(msa) :
    newMSA = {}
    for name in msa.keys() :
        newMSA[name] = []
        for seg in msa[name] :
            newMSA[name].append('')
    segLen = []
    for seg in msa.values()[0] :
        segLen.append(len(seg))
    nSegs = len(segLen)

    newSegLen = []
    for i in range(nSegs) :
        newLen = segLen[i]
        for col in range(segLen[i]) :
            nonEmptyCols = 0
            for name in msa.keys() :
                if msa[name][i][col] != '-' :
                    nonEmptyCols += 1
            if nonEmptyCols == 0:
                newLen -= 1
                continue
            for name in msa.keys() :
                newMSA[name][i] += msa[name][i][col]
        newSegLen.append(newLen)
    return newMSA

# Concatenate msa:
# From a map: seqName -> [[seg1], [seg2]]
# To a map: seqName -> whole seq, + [segment lengths]
def concatenate_msa(msa) : 
    segLen = []
    for seg in msa.values()[0] :
        segLen.append(len(seg))
    
    alignment = {}
    for k in msa :
        tmpAlign = ''
        for seg in msa[k] :
            tmpAlign += copy.deepcopy(seg)
        alignment[k] = tmpAlign
    return alignment, segLen

# Split MSA:
# reverse process of concatenate msa
def split_msa(MSA, segLengs) :
    msa = copy.deepcopy(MSA)
    alignment = {}
    for k in msa :
        alignment[k] = []
        for l in segLengs :
            alignment[k].append(copy.deepcopy(msa[k][0:l]))
            msa[k] = msa[k][l:]
    return alignment