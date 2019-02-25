#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
#
# concatAlign.py -> This code concatenates sequences from distinct files
#                   according to an identifier (id =species, accession
#                   number,etc) found on sequence headers.
#
# Usage: python concatAlign.py workingDirectory output.fasta
#
# Release date: 12/03/2018
# Last update: 12/03/2018
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import os
from sys import *
from Bio import SeqIO


dir = argv[1] # working directory
output = argv[2] # output file name
filenames = [inFile for inFile in os.listdir(dir) if inFile.endswith(".fasta")] # list of all file (.fasta) to be concatenated

dicSeq = {}
for inFile in filenames:
    # print inFile
    for entry in SeqIO.parse(open(dir+inFile),'fasta'):
        id, seq = entry.id, entry.seq
        id = id.split('_')[1] # change the separator as appropriate
        if id not in dicSeq.keys():
            dicSeq[id] = []
            dicSeq[id].append(str(seq))
        else:
            dicSeq[id].append(str(seq))

# saving the output
outFile = open(dir + output, "w")
for i,s in dicSeq.items():
    outFile.write(">" + i + "\n" + "".join(s) + "\n")
outFile.close()

print("\nSequence file created!\nPlease check it on " + dir)
