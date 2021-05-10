rm(list=ls())
setwd("/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/")
funcPath   <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/functions/"
resultPath <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/results/"


library(boot); library(tmvtnorm); #install.packages("tmvtnorm")
library(huge) #install.packages("huge")
library(igraph) #install.packages("igraph")
library(BDgraph) #install.packages("BDgraph")
library(SpiecEasi) 
#install.packages("reshape")
library(plyr);library(ggplot2);library(reshape);library(gridExtra)








load("QMPtree.RData")

x  <- RMP
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
		se <- spiec.easi(QMP, method = 'mb', sel.criterion = 'stars', 
		                 lambda.max = 0.9, nlambda = 100, lambda.min.ratio = 0.005/0.9, 
		                 pulsar.params = list(rep.num = 50, thresh = 0.2, subsample.ratio = 0.8 ) )

    save.image("result.RMPspieceasi.RData")



