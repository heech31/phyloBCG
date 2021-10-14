library(ape)
library(phytools)
library(phylobase)
library(phyloseq)
library(text2vec)
library(ppcor)
library(boot)
library(Matrix)
library(truncnorm); #Univariate truncated normal
library(huge);
library(igraph)
library(BDgraph)
library(plyr);
library(ggplot2);
library(reshape);
library(gridExtra)
library(doParallel)
library(SPRING)
library(SpiecEasi)
library(tmvtnorm)
library(mvtnorm)
library(truncdist)
library(doRNG)
library(corrplot)
#library(TreeTools)

# tmvtnorm::checkSymmetricPositiveDefinite examines
# positive definiteness of a matrix using determinant, which is incorrect
# use is.positive.definite
source("./functions/checkPD.R")
environment( checkSymmetricPositiveDefiniteCustom ) <- asNamespace("tmvtnorm")
assignInNamespace( "checkSymmetricPositiveDefinite", checkSymmetricPositiveDefiniteCustom, ns="tmvtnorm" )
