#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
#
#   pruneBeastTree.py -> This code removes clades from a tree given a list of
#                   taxa in an newline delimited auxiliary file.
#
# Usage: python workingDirectory treeFile auxFile
#
# Release date: 28/11/2017
# Last update: 14/12/2017
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


from sys import *
from Bio import Phylo
from io import StringIO

dir = argv[1]
treeFile = argv[2]
auxFile = open(dir + argv[3], "r").readlines() # taxa to be removed

treeHandle = open(dir + treeFile, 'r').readlines()

def beast2nwk(treeFile):
    # to get numeric names, and create a dict with their clade names
    dicNames = {}
    for line in treeFile:
        if line.startswith('\t\t  '):
            num, spp = line.strip().replace(',', '').split()
            dicNames[num] = spp

    # to get the line containing the tree data
    treedata = ''
    for line in treeFile:
        sep = 'tree TREE1 = '
        if line.startswith(sep):
            treedata = line.split(sep)[1].strip()

    handle = StringIO(treedata)
    tree = Phylo.read(handle, "newick")

    # to rename clade names
    for clade in tree.find_clades():
        if str(clade.name) in dicNames.keys():
            clade.name = dicNames[clade.name]

    return tree

nwkTree = beast2nwk(treeHandle)

# if in listPrune, these targets are removed
for tax in auxFile:
    nwkTree.prune(target=tax.strip())


# save tree
Phylo.write([nwkTree], dir+treeFile.split(".")[0]+"_pruned.tree", 'nexus')
print('\nFile save as \'' + treeFile.split(".")[0] + "_pruned.tree")


# to fix potential problems of formatting
treeLines = []
for line in open(dir+treeFile.split(".")[0]+"_pruned.tree").readlines():
    treeLines.append(line)

outFile = open(dir + treeFile.split(".")[0] + "_pruned.tree", "w")
for line in treeLines:
    if line.startswith(" Tree tree1"):
        line = line.replace("):0.00000[", ")[")
        line = line.replace(";", ":0.00000;")
        outFile.write(line)
    else:
        outFile.write(line)
outFile.close()