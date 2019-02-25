#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
#
# pairwiseModeller.py -> This code generates structural homology models of
#                       protein-protein interactions (PPI) between two proteins.
#                       Two sequence datasets of sequences homologous to those in the
#                       template are required. Proteins from each dataset are associated 
#                       in a pairwise manner to generate all possible PPI structures of 
#                       proteins present in both datasets.
#
# Usage: python pairwiseModeller.py workingDirectory viralSeqDatabase.fasta hostSeqDatabase PDBtemplate
#
# Release date: 24/05/2017
# Last update: 12/03/2018
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import os
import sys
import pylab
import modeller
import itertools
from sys import *
from Bio import SeqIO
from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb


dir = argv[1] # working directory
vProt = argv[2] # dataset of viral sequences homologous to the template structural sequence
hProt = argv[3] # same as above, for host sequences
template = argv[4] # pdb file of the protein-protein interaction

vhSystem = [vProt, hProt]
twoSets = []
dicSeqs = {}

# creates pairwise sets of virus-host interactors
def pairwisePPI(vhSystem):
    for domSet in vhSystem:
        protRep = []
        fasta_sequences = SeqIO.parse(open(domSet),'fasta')
        for fasta in fasta_sequences:
            id, seq = fasta.id, fasta.seq
            prot = id.split('|')[0]
            ssp = id.split('|')[2]
            dicSeqs[prot] = fasta.seq

            if prot not in protRep:
                protRep.append(prot + "-" + ssp)
        twoSets.append(protRep)
        fasta_sequences.close()
pairwisePPI(vhSystem)

lTuProts = list(itertools.product(twoSets[0], twoSets[1]))

# create a .ali file containing a pair of interacting proteins from the host and viral sequence database
def createDirAli(lTuProts):
    for pPair in lTuProts:
        folder = pPair[0].split("-")[0] + "-" + pPair[1].split("-")[0] + "_" + pPair[0].split("-")[1] + "-" + pPair[1].split("-")[1]
        ppi = folder.split("_")[0]
        if os.path.isfile(dir + folder + "/" + ppi + ".ali"):
            pass
        else:
            os.system("mkdir %s%s" % (dir, folder))
            outputFile = open(dir + folder + "/" + ppi + ".ali",'w')

            if pPair[0].split("-")[0] in dicSeqs:
                s1 = dicSeqs[pPair[0].split("-")[0]]
            if pPair[1].split("-")[0] in dicSeqs:
                s2 = dicSeqs[pPair[1].split("-")[0]]
            outputFile.write(">" + "P1;" + ppi + "\n")
            outputFile.write("sequence" + ":" + ppi + ":1:A:" + str(len(s1)+len(s2)) + ":B::::" + "\n")
            outputFile.write(str(s1) + "/" + str(s2) + "*\n")
createDirAli(lTuProts)


# s3_align2d.py routine: aligning query and template sequences ___________
def align2d():
    for pPair in lTuProts:
        folder = pPair[0].split("-")[0] + "-" + pPair[1].split("-")[0] + "_" + pPair[0].split("-")[1] + "-" + pPair[1].split("-")[1]
        ppi = folder.split("_")[0]
        if os.path.isfile(dir + "/" + folder + "/" + ppi + "_" + template.split(".")[0] + '.ali'):
            # print "File align2d exists."
            pass
        else:
            os.chdir(dir + folder)
            sys.stdout = open("align2d_" + ppi + ".log", 'w')
            env = environ()
            aln = alignment(env)
            mdl = model(env, file=dir + template, model_segment=('FIRST:A','LAST:B'))
            aln.append_model(mdl, align_codes=template.split(".")[0], atom_files=template)
            aln.append(file=ppi + '.ali', align_codes=ppi)
            aln.align2d()
            aln.write(file=ppi + "_" + template.split(".")[0] + '.ali', alignment_format='PIR')
            aln.write(file=ppi + "_" + template.split(".")[0] + '.pap', alignment_format='PAP')
align2d()


# s4_model-single.py routine: protein modelling __________________________
def modelSingle():
    for pPair in lTuProts:
        folder = pPair[0].split("-")[0] + "-" + pPair[1].split("-")[0] + "_" + pPair[0].split("-")[1] + "-" + pPair[1].split("-")[1]
        ppi = folder.split("_")[0]
        if os.path.isfile(dir + "/" + folder + "/" + ppi + '.ini'):
            pass
        else:
            os.chdir(dir + folder)
            os.system("cp ../%s ./" % (template))
            sys.stdout = open("model-single_" + ppi + ".log", 'w')
            env = environ()
            a = automodel(env, alnfile=ppi + "_" + template.split(".")[0] + '.ali',
                          knowns=template.split(".")[0], sequence=ppi,
                          assess_methods=(assess.DOPE,
                                          #soap_protein_od.Scorer(),
                                          assess.GA341))
            a.starting_model = 1
            a.ending_model = 5
            a.make()
modelSingle()


# This module searches for the best model based on their DOPE score
def findBestModel():
    for pPair in lTuProts:
        folder = pPair[0].split("-")[0] + "-" + pPair[1].split("-")[0] + "_" + pPair[0].split("-")[1] + "-" + pPair[1].split("-")[1]
        ppi = folder.split("_")[0]
        if os.path.isfile(dir + "/" + folder + "/" + "bestModel_"+ppi+".pdb"):
            pass
        else:
            os.chdir(dir + folder)
            minDope = 0
            bestPdb = ""
            listDopeSim = []

            for line in open("model-single_" + ppi + ".log").readlines()[-26:]:
                if line.startswith(ppi):
                    entry = line.split()[0] + "_" + line.split()[2]
                    listDopeSim.append(entry)
            for element in listDopeSim:
                if float(element.split("_")[1]) < minDope:
                    minDope = float(element.split("_")[1])
                    bestPdb = element.split("_")[0]
            os.system("cp %s ./bestModel_%s.pdb" % (bestPdb, ppi))
findBestModel()


# s5_evaluate_model.py routine ___________________________________________
def evaluate_model():
    for pPair in lTuProts:
        folder = pPair[0].split("-")[0] + "-" + pPair[1].split("-")[0] + "_" + pPair[0].split("-")[1] + "-" + pPair[1].split("-")[1]
        ppi = folder.split("_")[0]
    if os.path.isfile(dir + "/" + folder + "/" + ppi + '.profile'):
        pass
    else:
        os.chdir(dir + folder)
        sys.stdout = open("evaluate_model_" + ppi + ".log", 'w')
        model = 'bestModel_' + ppi + ".pdb"
        fileName = ppi + '.profile'

        log.verbose()    # request verbose output
        env = environ()
        env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
        env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

        # read model file
        mdl = complete_pdb(env, model)

        # Assess with DOPE:
        s = selection(mdl)   # all atom selection
        s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=fileName,
                      normalize_profile=True, smoothing_window=15)
evaluate_model()

