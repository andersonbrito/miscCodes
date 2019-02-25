#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
#
# degapper.py  ->  This code removes gaps giving a maximum threshold (%)
#                   of gaps per site in a multiple sequence alignment.
#
# Usage: python degapper.py workingDirectory gapPercentage
#
# Release date: 21/12/2017
# Last update: 21/12/2017
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from Bio import SeqIO
from sys import *
import os

dir = argv[1] # working directory
filenames = [inFile for inFile in os.listdir(dir) if inFile.endswith(".fasta")] # list of all file (.fasta) to be degapped
percGap = int(argv[2]) # if not greater the this cutoff, keep it


for inFile in filenames:
    # generates a list of headers and sequences to be transposed
    dicSeqs = {}
    fasta_sequences = SeqIO.parse(open(dir + inFile), 'fasta')
    for fasta in fasta_sequences:
        id, seq = fasta.description, fasta.seq
        seq = str(seq).upper()
        dicSeqs[id] = list(seq)

    # filter out position with more gaps than the percentage of allowed gaps (percGap variable)
    newSeqs = []
    for pos in list(map(list, zip(*list(dicSeqs.values())))): # transpose the MSA
        if not pos.count('-')/len(pos) > percGap:
            newSeqs.append(pos)

    # save the degapped file
    outFile = open(dir + inFile.split('.')[0] + '_degap.aln', "w")
    for num, pos in enumerate(list(map(list, zip(*list(newSeqs))))):
        entry = str(">" + list(dicSeqs.keys())[num] + "\n" + ''.join(pos) + "\n")
        outFile.write(entry)

    print('\nDegapped file saved: ' + inFile.split('.')[0] + '_degap.aln')
