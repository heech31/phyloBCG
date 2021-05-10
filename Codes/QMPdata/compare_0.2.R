rm(list=ls())
setwd("/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/")
funcPath   <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/functions/"
resultPath <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/results/"
figPath    <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/figures/"

## Sumarizing the results of PhyloBCG, SPRING, Spiec-Easiß
library(boot); library(tmvtnorm);
library(huge); library(SpiecEasi)#install.packages("huge")
library(igraph) #install.packages("igraph")
library(BDgraph) #install.packages("BDgraph")
library(plyr);library(ggplot2);library(reshape);library(gridExtra)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols <- gg_color_hue(20)

library(corpcor) #install.packages("corpcor")

## Load partial correlation matrices
# Tree
load(paste(resultPath,"pcor.tree_4chains.RData",sep=""))
pcor.tr <- parCor
# SPRING
load(paste(resultPath,"pcor.sp.RData",sep=""))
pcor.sp <- pcor.K
# SPIEC-EASI
load(paste(resultPath,"pcor.seRMP.RData",sep=""))
pcor.se <- pcor.S

rownames(pcor.sp) = colnames(pcor.sp) = rownames(pcor.tr)
rownames(pcor.se) = colnames(pcor.se) = rownames(pcor.tr)

#colnames(pcor.tr)


#------------------------------------------------------------------------------------------
#-------------------BCG vs SPRING ---------------------------------------------------------
#------------------------------------------------------------------------------------------
# Create a matrix to draw tile figure.
# SPRING is in lower matrix and BCG in upper matrix.
ltriFlag <- lower.tri(pcor.tr)
utriFlag <- upper.tri(pcor.tr)

TileMat <- matrix(0, ncol = dim(pcor.tr)[2], nrow = dim(pcor.tr)[1])

TileMat[ltriFlag] <- as.matrix(pcor.tr)[ltriFlag]

TileMat[utriFlag] <- as.matrix(pcor.sp)[utriFlag]

colnames(TileMat) <- colnames(pcor.tr); rownames(TileMat) <- colnames(pcor.tr)


# Create a coordinate vectors for geom_tile
# Get the locations of non-zero elements
##########################################################################################################
pcorthresh = 0.20 # threshold for absolute partial correlation
##########################################################################################################
isover_tr <- apply(pcor.tr, 1, function(x) any(abs(x) > pcorthresh))
isover_sp <- apply(pcor.sp, 1, function(x) any(abs(x) > pcorthresh))
isover_se <- apply(pcor.se, 1, function(x) any(abs(x) > pcorthresh))

isover_thr <- which( isover_tr | isover_sp |isover_se  )
length(isover_thr)


### This function returns partial correlations of three methods
compare.pcor = function(names){
  if(length(names)==1){
    corresponding.pcor <- cbind(pcor.tr[names[1],],pcor.sp[names[1],],pcor.se[names[1],])
    colnames(corresponding.pcor) <- c("Tree","SP","SE")
  }else{
    corresponding.pcor <- cbind(pcor.tr[names[1],names[2]],pcor.sp[names[1],names[2]],pcor.se[names[1],names[2]])
    colnames(corresponding.pcor) <- c("PhyloBCG","SPRING","SPIEC-EAIS")
    rownames(corresponding.pcor) <- paste(names[1],", ",names[2],sep="")
  }
  return(corresponding.pcor)
}




tab1 <- 
  rbind(
    compare.pcor(c("Dialister","Phascolarctobacterium")),
    compare.pcor(c("Oscillospira","Ruminococcus")),
    compare.pcor(c("Mitsuokella","Prevotella")),
    compare.pcor(c("Ruminococcus","Blautia")),
    ## PhyloBCG, SPRING
    compare.pcor(c("Oscillospira","Butyricimonas")),
    compare.pcor(c("Eubacterium","Peptococcus")),
    compare.pcor(c("Bacteroides","Bilophila")),
    compare.pcor(c("Akkermansia","Methanobrevibacter")),
    compare.pcor(c("Blautia","Methanobrevibacter")),
    compare.pcor(c("Prevotella","Bacteroides")),
    compare.pcor(c("Veillonella","Streptococcus")),
    compare.pcor(c("Bifidobacterium","Holdemania"))
  )

round(tab1,3)
