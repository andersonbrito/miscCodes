#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
#
# pruneNewickTree.py -> This code removes clades from a tree given a list of
#                   taxa in an auxiliary file.
#
# Usage: python directory auxFile listOfFiles
#
# Release date: 28/11/2017
# Last update: 14/12/2017
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


from Bio import Phylo
# from cStringIO import StringIO

dir = "/Users/anderson/boxsync/PhD/works/phylog/gene_trees/host/rhg7_n1234Phyml/"
treeFile = 'nect1234_phyml_rr.tree'
auxFile = open(dir + "listPrune.txt", "r").readlines() # taxa to remove

tree = Phylo.read(dir+treeFile, 'newick')

# to check clade names
for clade in tree.find_clades():
    if str(clade.name) != 'None':
        print(str(clade.name))


# list of taxa to prune from tree
listPrune = []
for tax in auxFile:
    listPrune.append(tax.strip())
    # print(tax)

# if in listPrune, these targets are removed
for c in listPrune:
    tree.prune(target=c)

print(tree)


# save tree
Phylo.write([tree], dir+treeFile.split(".")[0]+"_pruned.tree", 'newick')
print('\nFile save as \'' + treeFile.split(".")[0] + "_pruned.tree")


# to fix potential problems of formatting
treeLines = []
for line in open(dir+treeFile.split(".")[0]+"_pruned.tree").readlines():
    treeLines.append(line)

outFile = open(dir + treeFile.split(".")[0] + "_pruned.tree", "w")
for l in treeLines:
    line = line.replace("):", ")")
    # line = line.replace(")0.00000", "):0.3", 1000)
    # line = line.replace("0.00000", "0.3", 1000)
    outFile.write(line)
outFile.close()
