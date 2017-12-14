from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo

#Molluscan is already aligned
#Grab our fasta file containing our sequences
seq_record = SeqIO.parse("molluscan.fasta", "fasta")

"""
#For parsing through our sequences
for sequence in seq_record:
    print(sequence.id)
    print(repr(sequence.seq))
    print(len(sequence))
"""

#Used for if we need to access alignments as a list
#seq_record_list = list(seq_record)

#Perform msa here
#See 6.4 of biopython tutorial

#Does not work if we first parsed through sequence
#Output our aligned sequence as a phylip file
SeqIO.write(seq_record, "test.phy", "phylip")

#.phy seems to be main file type people work with when dealing with multi seq alignments
#Read our phy multiple alignment file
aln = AlignIO.read('test.phy', 'phylip')
print aln

#Using model "identity", scoring matrix to calculate distance
calculator = DistanceCalculator('identity')
#Grab distance matrix of a given alignment object
dm = calculator.get_distance(aln)
print(dm)

#Use neighbor joining to calculate our tree
constructor = DistanceTreeConstructor(calculator, 'nj')
tree = constructor.build_tree(aln)
print(tree)

#Write out our tree as a phyloxml file
Phylo.write(tree, 'testtree.xml', 'phyloxml')

tree2 = Phylo.read('testtree.xml',  'phyloxml')

#Draws trees
#Flip branches so deeper clades are displayed at top
tree.ladderize()
Phylo.draw(tree2)
