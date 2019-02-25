#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
#
# fetchGenomes.py -> This code performs bulk fetch of genome sequences from
#                   NCBI database. The auxiliary file contains a list of
#                   accession number separated by newline characters.
#
# Usage: python fetchGenomes.py workingDirectory auxiliaryFile
#
# Release date: 23/12/2017
# Last update: 12/03/2018
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from sys import *
from Bio import Entrez
try:
    import urllib.request as urllib2
except ImportError:
    import urllib2

Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are

dir = argv[1] # working directory
auxFile = open(dir + argv[2], "r").readlines() # list of accession numbers

for line in auxFile:
    accNo = line.strip()
    print(accNo)
    net_handle = Entrez.efetch(db="nucleotide", id=accNo, rettype="fasta", retmode="text")
    out_handle = open(dir + accNo + ".fasta", "w")
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()
    print("Saved " + accNo)
