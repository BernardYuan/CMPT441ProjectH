from Bio import Phylo
from cStringIO import StringIO

#Open the file and process every line
cfile = open('result.txt', 'r')
treeSet = []
for line in cfile.readlines():
    start = line.find('(')
    line = line[start:]
    while line.find('Inner') > 0:
        x = line.find('Inner')
        line = line[:x] + line[x+6:]
    treeSet.append(line)
print treeSet
for i in range(len(treeSet)):
    treeStr = treeSet[i]
    handle = StringIO(treeStr)
    tree = Phylo.read(handle, "newick")
    x = Phylo.draw(tree, branch_labels=lambda c: c.branch_length)
    raw_input()
       
cfile.close
