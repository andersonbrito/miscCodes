#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
#
# genomeStats.py -> This code searches for stats from viral genomes
#                   available on NCBI website, and retrieves these data
#                   in TSV format. Additional genetic information is
#                   appended at the end of the document, such as Average,
#                   Standard deviation, Minimum, Maximum genome size and
#                   number of proteins).
#
# Usage: python genomeStats.py taxid >> fileName.tsv
#
# Release date: 26/02/2015
# Last update: 02/03/2015
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from bs4 import BeautifulSoup as BS
import urllib2
import numpy as NP
import re
from sys import *

idNumber = argv[1]
outfile = open(str(idNumber) + '_genomeStats.txt', 'w')

webPage = urllib2.urlopen("http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=" + idNumber)
#
webPageData = webPage.read()
soup = BS(webPageData)
virTax = soup.find("i").text

dataTable = soup.find('table', class_='tblTxt')

if dataTable is None:
    outfile.write("There is no complete genome available for %s (taxid: %s).\n" % (virTax, idNumber))
    exit()
dataRows = dataTable.find_all('tr')[1:]

outfile.write("Taxonomic level: " + virTax + "\n")
outfile.write('\t'.join(("Taxid", "Species", "Segments",  "Genome size", "# proteins", '\n')))

lstgenSize = []
lstproNumber = []
lstsegNumber = []
allEntries = re.findall("\xc2\xa0nt</td>", str(soup))

for row in dataRows:
    rowCells = row.find_all('td')
    if (len(rowCells) < 10) or re.search('colspan', str(rowCells)):
        continue
    url = str(rowCells[0].find('a')['href'])
    rgx = re.search(r'id\=(\d+)', url)
    if rgx:
        accNo = rgx.group(1).strip()
    sppName = rowCells[0].text.strip()
    segNumber = rowCells[3].text

    if segNumber == "-":
        segNumber = "1"
    lstsegNumber.append(int(segNumber))
    if NP.std(lstsegNumber) > 0:
        segments = "Members of this taxa have %d-%d genomic segments" % (NP.min(lstsegNumber), NP.max(lstsegNumber))
    else:
        segments = "Members of this taxa have %d genomic segment(s)" % (NP.round(NP.mean(lstsegNumber)))
    genSize = rowCells[4].text.strip().encode('ascii', errors='ignore').replace("nt", "")
    lstgenSize.append(int(genSize))
    proNumber = rowCells[5].text.strip()
    if proNumber != "-":
        lstproNumber.append(int(proNumber))

        outfile.write('\t'.join((accNo, sppName, segNumber, genSize, proNumber, '\n')))


outfile.write("\nStats from %d complete genomes of %s.\n%s:\n" % (len(lstgenSize), virTax, segments))
outfile.write('\t'.join(("Data", "Avg", "Std", "Min", "Max", '\n')))
outfile.write('\t'.join(("Genome", str(NP.round(NP.mean(lstgenSize))), str(NP.round(NP.std(lstgenSize))), str(NP.min(lstgenSize)), str(NP.max(lstgenSize)), '\n')))
outfile.write('\t'.join(("Proteins", str(NP.round(NP.mean(lstproNumber))), str(NP.round(NP.std(lstproNumber))), str(NP.min(lstproNumber)), str(NP.max(lstproNumber)), '\n\n')))

if len(lstgenSize) != len(allEntries):
    absGenones = len(allEntries) - len(lstgenSize)
    print("A total of %d were found, however, NCBI has %d complete genome(s) of ""%s"". \n" \
          "Please access the following website to check the %d absent genomes: http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=%s\n" % (len(lstgenSize), len(allEntries), virTax, absGenones, idNumber))

print('Done!\nResults saved in your working directory: ' + str(idNumber) + '_genomeStats.txt')
