library(devtools)
install_github("mattheusmarzola/GOpruning")
#source("https://bioconductor.org/biocLite.R")

library(org.EcK12.eg.db)
library(STRINGdb)
library(igraph)
library(RedeR)
library(GO.db)
library(plyr)
library(data.tree)
library(RColorBrewer)
library(GOpruning)

initialTerm<- "GO:0051090"
goHierarchy <- getHierarchyByTerm(rootTerm = initialTerm, ontology = 'BP')



