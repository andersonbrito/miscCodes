# Script for creating plot with multiple scatterplots
# Written by: Anderson Brito
# Date: 19 Feb 2016
# Contact: andersonfbrito@gmail.com
# R version: R version 3.3.2

setwd("/path/to/working/directory/")

dope_sim <- read.csv("example_multiscatter.tsv", header=TRUE, sep="\t")

dope_sim

virIntPos <- dope_sim$virIntPos
virDomPos <- dope_sim$virDomPos
virIdent <- dope_sim$virIdent
hostIntPos <- dope_sim$hostIntPos
hostDomPos <- dope_sim$hostDomPos
hostIdent <- dope_sim$hostIdent
dope_score <- dope_sim$dope_score
rms <- dope_sim$rms
infection <- as.numeric(dope_sim$infection)

d <- data.frame(virDomPos=virDomPos, hostDomPos=hostDomPos, dope_score=dope_score, rms=rms)
pairs(d)