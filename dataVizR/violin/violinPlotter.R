# Script for creating violin plots
# Written by: Anderson Brito
# Date: 03 May 2016
# Contact: andersonfbrito@gmail.com
# R version: R version 3.3.2

library(ggplot2)

setwd("/path/to/working/directory/")

erelgv <- read.csv("example_violin.txt", header=TRUE, sep="\t")

vir <- as.factor(erelgv$vir)
div <- erelgv$div
cat <- erelgv$cat

# Basic violin plot
p <- ggplot(erelgv, aes(x=vir, y=div), size=2) + 
  geom_violin(trim=FALSE)  +
  scale_y_continuous(breaks=seq(0, 12, 1)) +
  labs(x="Virus", y = "Number of nonsynonymous substitutions") +
  geom_jitter(shape=16, position=position_jitter(0.25), aes(colour = factor(cat)),
              size=I(2.5))

# # Rotate the violin plot
# p + coord_flip()

p  + theme_bw() + theme(legend.position="none") +
  scale_colour_manual(values = c("#bbddffff", "#5d90d4ff", "#0044aaff", "#6bb200ff", "#006600ff", "#ffcc00ff"))
