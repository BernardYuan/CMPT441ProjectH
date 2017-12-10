Currently the code takes in a FASTA file, assumes it is already aligned, computes the 
distance matrix and creates a tree using NJ, then  outputs the tree ast testtree.xml


Required modules
-numpy
-scipy
-matplotlib
-biopython

Phylo handles distance matrix, base phylo tree construction from NJ

-Guide for IO of FASTA and PHY
http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc47

-Guide for matrix and tree
http://biopython.org/wiki/Phylo

If MSA needed from FASTA, use MUSCLE

File types
Sequences.fasta
AlignedSequences.phy
Trees.xml