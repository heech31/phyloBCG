rm(list=ls())
setwd("/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/")
funcPath   <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/functions/"
resultPath <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/results/"


library(boot); library(tmvtnorm); #install.packages("tmvtnorm")
library(huge) #install.packages("huge")
library(igraph) #install.packages("igraph")
library(BDgraph) #install.packages("BDgraph")
library(SPRING) 
#install.packages("reshape")
library(plyr);library(ggplot2);library(reshape);library(gridExtra)








load("QMPtree.RData")

x  <- QMP
np <- dim(x)
n  <- np[1] # Sample size
p  <- np[2]  # Data dimension
K <- 2   # Latent space dimension


colnames(x)
dim(x)
qr(x)$rank






########################################################################
########   Parameters for precision matrix sampling ####################
########################################################################


		bt <- Sys.time()
		set.seed( 77843 )		
		sp <- SPRING(QMP, Rmethod = "approx", quantitative = TRUE, thresh = 0.2, verbose = FALSE,
		             lambdaseq = "data-specific", nlambda = 100, rep.num = 50)
		
    save.image("result.QMPsping.RData")



