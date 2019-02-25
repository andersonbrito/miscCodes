# Script for creating maps with pinpoints given a list of coordinates
# Written: Anderson Brito
# Date: 15 March 2016
# Contact: andersonfbrito@gmail.com
# R version: R version 3.3.2

library(maps)

setwd("/Users/anderson/boxsync/PhD/courses/math/misc/maps/")

# file containing GPS points in decimal degrees
samps <- read.csv("coord.csv")

# creating the map
map('worldHires', 'Brazil', col="gray80", fill=TRUE)
map('worldHires', '', col="gray95", fill=TRUE, add=TRUE)

# plotting the pinpoints
points(samps$lon, samps$lat, pch=19, col="firebrick", cex=1)
