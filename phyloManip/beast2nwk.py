#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Python version: Python 3
#
# concatAlign.py -> This Script converts Beast trees into Newick format.
#
# Usage: python beast2nwk.py workingDirectory beastTree
#
# Release date: 14/12/2017
# Last update: 12/03/2018
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from sys import *
from Bio import Phylo
from io import StringIO
import re

dir = argv[1]
inFile = argv[2]
treeFile = open(dir + inFile, 'r').readlines()

# to get numeric names, and create a dict with their clade names
dicNames = {}
for line in treeFile:
    if line.startswith('\t\t  '):
        num, spp = line.strip().replace(',', '').split()
        dicNames[num] = spp

print(dicNames)

# to get the line containing the tree data
treedata = ''
for line in treeFile:
    sep = 'treetree1='
    if sep.lower() in line.lower().replace(' ',''):
        print(line)
        line = line.replace(' ','').replace('Tree', 'tree').replace('TREE', 'tree')
        treedata = line.replace(' ','').split(sep)[1].strip()

handle = StringIO(treedata)
tree = Phylo.read(handle, "newick")
# print(tree)


# to rename clade names
for clade in tree.find_clades():
    if str(clade.name) in dicNames.keys():
        clade.name = dicNames[clade.name]

    if str(clade.name) == 'None':
        listComm = str(clade.comment).split(",")
        for c in listComm:
            if 'posterior' in c:
                posterior = c.split("=")[-1]
                clade.confidence = float(posterior) * 100
        clade.comment = ""
    if str(clade.name) != 'None':
        clade.comment = ""

print(tree)

Phylo.write([tree], dir + inFile.split(".")[0] + "_nwk.tree", 'newick')
print('\nFile save as \'' + inFile.split(".")[0] + "_nwk.tree'")


# to fix potential problems of formatting
treeLines = []
for line in open(dir + inFile.split(".")[0] + "_conv.tree").readlines():
    treeLines.append(line)

outFile = open(dir + inFile.split(".")[0] + "_conv.tree", "w")
for l in treeLines:
    line = line.replace("):", ")")
    outFile.write(line)
outFile.close()
