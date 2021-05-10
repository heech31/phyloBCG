rm(list=ls())
setwd("/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/")
funcPath   <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/functions/"
resultPath <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/results/"
figPath    <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/figures/"


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
#load(paste(resultPath,"pcor.tree.RData",sep=""))
#load(paste(resultPath,"pcor.tree_1e-3.RData",sep=""))
#load(paste(resultPath,"pcor.tree_z.RData",sep=""))
load(paste(resultPath,"pcor.tree_4chains.RData",sep=""))
#load(paste(resultPath,"pcor.tree.RData",sep=""))
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


# create a coordinate vectors for geom_tile
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


colnames(pcor.tr)[isover_tr]
colnames(pcor.se)[isover_se]

compare.pcor(c("Dorea","Collinsella"))
compare.pcor(c("Blautia","Succinivibrio"))
compare.pcor(c("Bifidobacterium","Dorea"))
compare.pcor(c("Parabacteroides","Holdemania"))
compare.pcor(c("Peptococcus","Megamonas"))


compare.pcor(c("Lactobacillus","Enterococcus"))
compare.pcor(c("Lactobacillus","Streptococcus"))
compare.pcor(c("Lactobacillus","Haemophilus"))








compare.pcor(c("Parabacteroides","Holdemania"))
compare.pcor(c("Akkermansia","Methanobrevibacter"))

plot(graph_from_adjacency_matrix(tmp.tr*1,"undirected"))
#plot(graph_from_adjacency_matrix(tmp.se*1,"undirected"))
#plot(graph_from_adjacency_matrix(tmp.sp*1,"undirected"))



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

compare.pcor(c("Roseburia","Paraprevotella"))

which( colnames(x) == "Dialister" )
which( colnames(x) == "Phascolarctobacterium" )
plot(x[,c("Dialister","Phascolarctobacterium")])
plot(log(x[,c("Dialister","Phascolarctobacterium")]+1))
plot(zhat[,c(33,39)])


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
    
    
    compare.pcor(c("Prevotella","Bacteroides")),
    compare.pcor(c("Anaerostipes","Methanobrevibacter")),
    compare.pcor(c("Veillonella","Streptococcus")),  
    compare.pcor(c("Bifidobacterium","Holdemania"))
  )

round(tab1,3)






compare.pcor(c("Ruminococcus","Blautia"))
compare.pcor(c("Eubacterium","Peptococcus"))

compare.pcor(c("Oscillospira","Peptococcus"))
compare.pcor(c("Oscillospira","Eubacterium"))


compare.pcor(c("Butyricimonas","Peptococcus"))
compare.pcor(c("Butyricimonas","Eubacterium"))







# A helps M's growth
# Bac release lactate and A uses lactate to produce butyrate

#compare.pcor(c("Akkermansia","Bacteroides"))

###############################################################
### ALL
compare.pcor(c("Dialister","Phascolarctobacterium"))
compare.pcor(c("Oscillospira","Ruminococcus"))
compare.pcor(c("Haemophilus","Veillonella"))
## PhyloBCG, SPRING
compare.pcor(c("Mitsuokella","Prevotella"))
compare.pcor(c("Ruminococcus","Blautia"))
compare.pcor(c("Akkermansia","Methanobrevibacter"))
compare.pcor(c("Blautia","Methanobrevibacter"))

compare.pcor(c("Eubacterium","Peptococcus"))
compare.pcor(c("Bilophila","Bacteroides"))
## PhyloBCG, SPIEC-EASI
compare.pcor(c("Blautia","Methanobrevibacter"))
compare.pcor(c("Prevotella","Bacteroides"))


## PhyloBCG
compare.pcor(c("Streptococcus","Veillonella"))
compare.pcor(c("Streptococcus","Dialister"))
compare.pcor(c("Veillonella","Dialister"))

compare.pcor(c("Christensenella","Oscillospira"))








##########################################################################################################
pcorthresh = 0.1 # threshold for absolute partial correlation
##########################################################################################################
isover_tr <- apply(pcor.tr, 1, function(x) any(abs(x) > pcorthresh))
compare.pcor(c(colnames(pcor.tr)[isover_sp][7]))
colnames(pcor.tr)[isover_sp][7]


pp=0
pp= pp +1
compare.pcor(c(   colnames(pcor.tr)[isover_tr][pp]  ))
colnames(pcor.tr)[isover_tr][pp]

# Bacteroides and Sutterella
# Butyricimonas and Sutterella
# Bifidobacterium vs Coprobacillus (+), Holdemania (-), Turicibacter (+)
# Actinomyces vs Sutterella (-)
# Collinsella vs Adlercreutzia (+)
# Aldercreutzia vs Catenibacterium (-)
# Slackia vs Coprobacillus (-)
# Anaerostipes vs Dehalobacterium (+)
# Akkermansia vs Ruminococcus (+)


colnames(pcor.tr)[isover_sp]


compare.pcor(c("Actinomyces"))
compare.pcor(c("Bacteroides"))
compare.pcor(c("Ruminococcus"))
compare.pcor(c("Acidaminococcus"))


compare.pcor(c("Christensenella","Methanobrevibacter"))
compare.pcor(c("Haemophilus"))
compare.pcor(c("Streptococcus"))
compare.pcor(c("Enterococcus"))
compare.pcor(c("Peptococcus"))



compare.pcor(c("Oscillospira","Dehalobacterium"))

## SE
tab2 <- 
rbind(
compare.pcor(c("Roseburia","Anaerostipes")),
compare.pcor(c("Ruminococcus","Faecalibacterium")),
compare.pcor(c("Ruminococcus","Oscillospira")),
compare.pcor(c("Coprococcus","Blautia")),
compare.pcor(c("Roseburia","Faecalibacterium")),
compare.pcor(c("Roseburia","Faecalibacterium")),
compare.pcor(c("Sutterella","Faecalibacterium"))
)

round(tab2,3)

compare.pcor(c("Coprococcus"))
compare.pcor(c("Faecalibacterium"))

compare.pcor(c("Mitsuokella","Prevotella"))
compare.pcor(c("Ruminococcus","Blautia"))

compare.pcor(c("Ruminococcus","Prevotella"))
compare.pcor(c("Mitsuokella","Blautia"))
compare.pcor(c("Mitsuokella","Ruminococcus"))


## SPRING
tab3 <- 
rbind(
compare.pcor(c("Parabacteroides","Bacteroides")),
compare.pcor(c("Butyricimonas","Oscillospira")),
compare.pcor(c("Butyricimonas","Odoribacter")),
compare.pcor(c("Dorea","Blautia"))
)
round(tab3,3)



compare.pcor(c("Akkermansia"))
compare.pcor(c("Blautia","Methanobrevibacter"))


c("Roseburia",NULL)








coord <- cbind(rowv = rep(rownames(pcor.tr)[isover_thr], each = length(isover_thr)), 
               colv = rep(rownames(pcor.tr)[isover_thr], length(isover_thr)))

df_Tilemat <- cbind.data.frame(coord, partialcorr = c(TileMat[isover_thr,isover_thr]))

# rows start from bottom to top and cols start from left to right.
# but we want start from the left top corner for both triangles. So make colv factor in reverse level.
df_Tilemat$rowv <- factor(df_Tilemat$rowv, levels = rownames(pcor.tr)[isover_thr])
df_Tilemat$colv <- factor(df_Tilemat$colv, levels = rev(rownames(pcor.tr)[isover_thr]))



p02 <- length(isover_thr)
tile1 <- ggplot(df_Tilemat, aes(x = rowv, y = colv), fill = partialcorr)        ## global aes

tile1 <- tile1 + geom_tile(aes(fill = partialcorr), color = "black", alpha = 1, width = 1, height = 1) ## to get the rect filled
tile1 <- tile1 + scale_fill_gradient2(high = "blue", mid = "white", low = "red", name = "Partial Correlation", limit=c(-0.5,0.5))## color of the corresponding aes
tile1 <- tile1 + coord_cartesian(xlim = c(1, p02), ylim = c(1, p02), clip = 'off') 
tile1 <- tile1 + theme_bw()

tile1 <- tile1 +  
        theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0), 
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.position = c(-0.05, 1.25), legend.justification="left", legend.direction="horizontal", plot.margin = margin(1, 1, 1, 0, "cm")) 

tile1 <- tile1 + annotate("text", x = (p02+2)/2, y = -0.7, label = "BCG on Quantitative counts") +
        annotate("text", x = p02+1.5, y = (p02+2)/2, angle = 90, label = "SPRING on Quantitative counts") 

tile1 <- tile1 +  scale_x_discrete(position = "top") + 
                  scale_alpha_continuous(guide=FALSE) +  scale_size_continuous(guide=FALSE) + 
                  geom_segment(aes(x=0.5, y=(p02+0.5), xend = (p02+0.5), yend=0.5), color="grey", size = 0.5)
tile1

pdf(file = paste(figPath,"Tile_BCG_SPRING_20.pdf",sep=""), width = 8, height = 8)
tile1
dev.off()











#------------------------------------------------------------------------------------------
#-------------------BCG vs SPIEC-EASI  ----------------------------------------------------
#------------------------------------------------------------------------------------------
# Create a matrix to draw tile figure.
# SPRING is in lower matrix and BCG in upper matrix.
ltriFlag <- lower.tri(pcor.tr)
utriFlag <- upper.tri(pcor.tr)

TileMat <- matrix(0, ncol = dim(pcor.tr)[2], nrow = dim(pcor.tr)[1])

TileMat[ltriFlag] <- as.matrix(pcor.tr)[ltriFlag]

TileMat[utriFlag] <- as.matrix(pcor.se)[utriFlag]

colnames(TileMat) <- colnames(pcor.tr); rownames(TileMat) <- colnames(pcor.tr)


# create a coordinate vectors for geom_tile
# Get the locations of non-zero elements



isover_tr <- apply(pcor.tr, 1, function(x) any(abs(x) > pcorthresh))
isover_sp <- apply(pcor.sp, 1, function(x) any(abs(x) > pcorthresh))
isover_se <- apply(pcor.se, 1, function(x) any(abs(x) > pcorthresh))

isover_thr <- which( isover_tr | isover_sp |isover_se  )
# isover_thr <- 1:dim(pcor.sp)[2]
coord <- cbind(rowv = rep(rownames(pcor.tr)[isover_thr], each = length(isover_thr)), 
               colv = rep(rownames(pcor.tr)[isover_thr], length(isover_thr)))

df_Tilemat <- cbind.data.frame(coord, partialcorr = c(TileMat[isover_thr,isover_thr]))

# rows start from bottom to top and cols start from left to right.
# but we want start from the left top corner for both triangles. So make colv factor in reverse level.
df_Tilemat$rowv <- factor(df_Tilemat$rowv, levels = rownames(pcor.tr)[isover_thr])
df_Tilemat$colv <- factor(df_Tilemat$colv, levels = rev(rownames(pcor.tr)[isover_thr]))



p02 <- length(isover_thr)
tile2 <- ggplot(df_Tilemat, aes(x = rowv, y = colv), fill = partialcorr)        ## global aes

tile2 <- tile2 + geom_tile(aes(fill = partialcorr), color = "black", alpha = 1, width = 1, height = 1) ## to get the rect filled
tile2 <- tile2 + scale_fill_gradient2(high = "blue", mid = "white", low = "red", name = "Partial Correlation", limit=c(-0.5,0.5))## color of the corresponding aes
tile2 <- tile2 + coord_cartesian(xlim = c(1, p02), ylim = c(1, p02), clip = 'off') 
tile2 <- tile2 + theme_bw()

tile2 <- tile2 +  
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0), 
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.position = c(-0.05, 1.25), legend.justification="left", legend.direction="horizontal", plot.margin = margin(1, 1, 1, 0, "cm")) 

tile2 <- tile2 + annotate("text", x = (p02+2)/2, y = -0.7, label = "BCG on Quantitative counts") +
  annotate("text", x = p02+1.5, y = (p02+2)/2, angle = 90, label = "SPIEC-EASI on Quantitative counts") 

tile2 <- tile2 +  scale_x_discrete(position = "top") + 
  scale_alpha_continuous(guide=FALSE) +  scale_size_continuous(guide=FALSE) + 
  geom_segment(aes(x=0.5, y=(p02+0.5), xend = (p02+0.5), yend=0.5), color="grey", size = 0.5)
tile2

pdf(file = paste(figPath,"Tile_BCG_SPIECEASI_20RMP.pdf",sep=""), width = 8, height = 8)
tile2
dev.off()



pdf(file = paste(figPath,"Tile_BCG_SRPING_SPIECEASI_20.pdf",sep=""), width = 16, height = 8)
grid.arrange(tile1,tile2,ncol=2)
dev.off()


pdf(file = paste(figPath,"Tile_BCG_SRPING_SPIECEASI_20_vertical.pdf",sep=""), width = 8, height = 16)
grid.arrange(tile1,tile2,ncol=1)
dev.off()






dim(TileMat[isover_thr,isover_thr])





