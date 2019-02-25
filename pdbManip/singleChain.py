#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
#
# singleChain.py -> This code splits a multichain PDB file into its
#                   multiple individual chains, saving them as output.
#
# Usage: python singleChain.py workingDirectory pdbFile
#
# Release date: 30/12/2017
# Last update: 30/12/2017
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


from Bio.PDB import PDBParser, PDBIO
from sys import *
import os

dir = argv[1]
inFile = argv[2]
newname = '_chain'

io = PDBIO()
pdb = PDBParser().get_structure(newname, dir + inFile)

for chain in pdb.get_chains():
    io.set_structure(chain)
    io.save(dir + inFile.split('.')[0] + newname + chain.get_id() + ".pdb")
