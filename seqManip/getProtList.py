#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
#
# getProtList.py -> This code downloads a list of proteins given a list
#                   of uniprot identifiers. The auxiliary file must
#                   contain tab separated lines with a uniprot id followed
#                   by a pfam code (ex. PF00136), or some other information.
#                   The pfam code is relevant only when running the script
#                   on 'domain' mode. In 'protein' mode, NCBI accession codes
#                   are also valid (to fetch from Entrez)
#
# Usage: python getProtList.py workingDirectory auxiliaryFile.tsv mode_proteinORdomain
#
# Release date: 19/05/2017
# Last update: 21/12/2017
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


try:
    import urllib.request as urllib2
except ImportError:
    import urllib2
from bs4 import BeautifulSoup as BS
import os
from Bio import SeqIO
from Bio import Entrez
from sys import *
Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are


dir = argv[1]
auxFile = argv[2]

protLst = []
for accno in open(dir+auxFile).readlines():
    protLst.append(accno)

# generate list of proteins from Uniref
mode = argv[3] # select 'protein' to fetch full sequence, or 'domain' to fetch only domain sequence (pfam code must be specified, see description above)
aaBefore = 0 # used on 'domain' mode only: if you want to include upstream flaking amino acids of the domain, define how long should this 'sequence envelope' be
aaAfter = 0 # same as above, but for amino acids downstream to the domain end position

# if the list of proteins are large, this module enables resuming a broken run
partList = []
outputFile = ""
if os.path.isfile(dir + auxFile.split(".")[0] + '.fasta'):
    fasta_sequences = SeqIO.parse(open(dir + auxFile.split(".")[0] + '.fasta'),'fasta')
    for fasta in fasta_sequences:
        id, seq = fasta.id, fasta.seq
        prot = id.strip().split('|')[1]
        if prot not in partList:
            partList.append(prot)
    fasta_sequences.close()
    outputFile = open(dir + auxFile.split(".")[0] + '.fasta','a')
else:
    outputFile = open(dir + auxFile.split(".")[0] + '.fasta','w')


# fetching domain or protein sequences from Uniprot or NCBI
allProt = []
c = 0
for item in protLst:
    item = item.strip()
    prot, dom = item.split('\t') # the auxiliary file must contain tab separated lines

    if prot in partList:
        pass
    else:
        try:
            link = urllib2.urlopen("http://www.uniprot.org/uniprot/{0}.fasta".format(prot.rstrip()))
            fastaSeq = BS(link,  "html.parser").string
            if not fastaSeq:
                pass
            else:
                seq = "".join(fastaSeq.split('\n')[1:])
        except:
            print(prot + ' not found!')
            handle = Entrez.efetch(db="protein", id=prot, rettype="gb")
            for seq_record in SeqIO.parse(handle, "gb"):
                seq = str(seq_record.seq)

        if mode is "protein":
            if prot not in allProt:
                print(prot)
                c += 1
                # print seq
                outputFile.write(">" + dom + "|" + prot + "\n")
                outputFile.write(seq + "\n")
                allProt.append(prot)
            else:
                print('duplicate →', prot)
                c += 1
                outputFile.write(">" + dom + "|" + prot + "\n")
                outputFile.write(seq + "\n")
                allProt.append(prot)

        if mode is "domain":
            pPage = urllib2.urlopen("http://pfam.xfam.org/protein/" + prot)
            pPageData = pPage.read()
            soup = BS(pPageData, "html.parser")
            pDataTable = soup.find('table', class_='resultTable details')
            if pDataTable:
                dataRows = pDataTable.find_all("tr")
                for row in dataRows:
                    if dom in str(row):
                        outputFile.write(">" + prot + "|" + dom + "\n")
                        if int(row.find_all("td")[2].text) - aaBefore < 0:
                            start = 1
                        else:
                            start = int(row.find_all("td")[2].text) - aaBefore
                        if int(row.find_all("td")[3].text) + aaAfter > len(seq):
                            end = len(seq)
                        else:
                            end = int(row.find_all("td")[3].text) + aaAfter
                        outputFile.write(str(seq[start:end] + "\n"))

print('Done!')
print(str(c), 'proteins fetched!')
outputFile.write("### END\n")
