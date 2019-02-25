#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
#
# repertsoup.py -> Provided an accession number of genomes (e.g. EU082088) or
#                   proteomes (e.g. UP000009136), this code retrieves domain
#                   repertoire from their proteins available on Pfam, and also
#                   gets the number of protein structures available for each
#                   domain from PDB.
#
# Usage: repertsoup.py 'workingDirectory' 'RefSeq or Uniprot Proteome Accno'
#
# Release date: 26/02/2015
# Last update: 02/03/2018
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from bs4 import BeautifulSoup as BS
try:
    import urllib.request as urllib2
except ImportError:
    import urllib2
import re
from sys import *
import os

dir = argv[1]
sAccNo = argv[2]

partList = []
listOutputs = []

outputFile = ""
if os.path.isfile(dir+sAccNo + '.txt'):
    for line in open(dir+sAccNo + '.txt').readlines()[1:]:
        rgx = re.search(r'(\w+)\t\w+\t.+', line)
        if rgx:
            prot = rgx.group(1).strip()
            if prot not in partList:
                partList.append(prot)
        listOutputs.append(line.strip())
    outputFile = open(dir+sAccNo + '.txt','a')
else:
    outputFile = open(dir+sAccNo + '.txt','w')


if sAccNo.startswith("UP"):
    link = "http://www.uniprot.org/uniprot/?query=proteome:{0}&format=list".format(sAccNo.rstrip())
else:
    link = "http://www.uniprot.org/uniprot/?query={0}&format=list".format(sAccNo.rstrip())

lPage = urllib2.urlopen(link)
listEdges = []


# This module searches for PDB entries of domains on the Pfam website, and returns the number of available structures.
def searchPDB(dAccNo):
    domPage = urllib2.urlopen("http://pfam.xfam.org/family/" + str(dAccNo))
    domPageData = domPage.read()
    soup = BS(domPageData, "html.parser")
    pdbEntries = str(soup.find('div', id="icons").find('span', id="structIcon").find('em').string).strip()

    return pdbEntries

# This module retrieve domain information for each protein in the genome or proteome
def searchDom(protList):
    speciesName = ""
    protDoms = []
    for pAccNo in protList:
        pPage = urllib2.urlopen("http://pfam.xfam.org/protein/" + pAccNo)
        pPageData = pPage.read()
        soup = BS(pPageData, "html.parser")
        pDataTable = soup.find('table', class_='resultTable details')

        if pDataTable is None:
            pPage = urllib2.urlopen("http://www.uniprot.org/uniprot/" + pAccNo)
            pPageData = pPage.read()
            soup = BS(pPageData, "html.parser")
            pDataTable = soup.find('table', class_='databaseTable DOMAIN')
            if pDataTable is None:
                pass
            else:
                if speciesName is "" and os.stat(dir+sAccNo + '.txt').st_size == 0:
                    speciesName = soup.find('div', id='content-organism').string
                    header = str("Results for " + sAccNo + " - " + speciesName)
                    outputFile.write(header + "\n")
                    print(header)

                dataRows = pDataTable.find_all('tr')
                if (len(dataRows) < 2) or not re.search('Pfam', str(dataRows)):
                    t = (pAccNo,'NA')
                    listEdges.append(t)
                    output = str(pAccNo + "\t" + 'Unknown' + "\t" + '-')
                    if output not in listOutputs:
                        listOutputs.append(output)
                        outputFile.write(output + "\n")
                        print(output)
                        continue

                for row in dataRows:
                    if (re.search('Pfam', str(row))):
                        dataCols = row.find_all(href=re.compile("family"))
                        for item in dataCols:
                            url = item.get('href')
                            rgx = re.search(r'family/(\w+)', url)
                            if rgx:
                                dAccNo = rgx.group(1).strip()
                                t = (pAccNo,str(dAccNo))
                                listEdges.append(t)
                                output = str(pAccNo + "\t" + dAccNo + "\t" + searchPDB(dAccNo))
                                if output not in listOutputs:
                                    listOutputs.append(output)
                                    outputFile.write(output + "\n")
                                    print (output)
        else:
            if speciesName is "" and os.stat(dir+sAccNo + '.txt').st_size == 0:
                speciesName = soup.find('tbody').find_all('a')[0].text.strip()
                header = str("Results for " + sAccNo + " - " + speciesName)
                outputFile.write(header + "\n")
                print(header)

            dataRows = pDataTable.find_all("td")
            if not (re.search('Pfam', str(dataRows))):
                output = str(pAccNo + "\t" + 'Unknown' + "\t" + '-')
                if output not in listOutputs:
                    listOutputs.append(output)
                    outputFile.write(output + "\n")
                    print(output)
                continue
            for row in dataRows:
                if (re.search('Pfam', row.text)):
                    dAccNo = row.get('class')[0].split('_')[1]
                    output = str(pAccNo + "\t" + dAccNo + "\t" + searchPDB(dAccNo))
                    if output not in listOutputs:
                        listOutputs.append(output)
                        outputFile.write(output + "\n")
                        print(output)
    return protDoms

alterLink = ""

fullList = []
if alterLink == "":
    for pAccNo in lPage.read().splitlines():
        pAccNo = str(pAccNo.decode('utf-8'))
        fullList.append(pAccNo)
        fullList = sorted(fullList)
else:
    lPage = urllib2.urlopen(alterLink)
    for pAccNo in lPage.read().splitlines():
        pAccNo = str(pAccNo.decode('utf-8'))
        fullList.append(pAccNo)
        fullList = sorted(fullList)

try:
    protList = [prot for prot in fullList if prot not in partList]
    protList.insert(0, partList[-1])
    protList = sorted(protList)
except:
    pass

searchDom(protList)

outputFile.write("END\n")
outputFile.close()