rm(list=ls())
setwd("/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/")
funcPath   <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/functions/"
resultsPath <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/ada/adaQMP/"
figPath    <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/figures/"

library(boot); library(tmvtnorm); library(coda)
library(huge) #install.packages("huge")
library(igraph) #install.packages("igraph")
library(BDgraph) #install.packages("BDgraph")
library(plyr);library(ggplot2);library(reshape);library(gridExtra)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols <- gg_color_hue(20)

library(corpcor) #install.packages("corpcor")
load('~/Dropbox/Research/TAMU/graphical/QMPdata/QMPtree.RData')
x <- QMP


likelihood <- function(z,R){
  np <- dim(z)
  n <- np[1]
  p <- np[2]
  log.likes <- rep(0,n)
  for(ii in 1:n){
    log.likes[ii] <- dmvnorm(z[ii,],mean=rep(0,p), sigma=R, log=TRUE)
  }
  log.lik <- sum(log.likes)
  list(log.lik=log.lik, log.likes=log.likes)
}


##################################################################################
##################################################################################
##################################################################################
## Thin sample by 5, 30,000/nthin
nmc <- 10000
p   <- 54
nchain <- 4
nthin   <- 4
thinInd <-  nthin*(1:(nmc/nthin) )
ngibbs  <- length( nthin*(1:(nmc/nthin)))
lp5 <- matrix(0,nchain,ngibbs)
pi_gibbs <- array(0,dim=c(p,p,ngibbs*nchain))
C_gibbs  <- array(0,dim=c(p,p,ngibbs*nchain))
sig2_gibbs <- rep(0,ngibbs*nchain)

dim(gibbsSample$U_gibbs)


# Posterior latent mean position
T_post_mean <- vector("list",4)



#ch=1
for(ch in 1:nchain){
  gibbsSample <- readRDS( paste(resultsPath,paste("result.QMPchain_",ch,".rds",sep=""), sep="" ) )

  print(ch)
  lp <- rep(0,ngibbs)
  for(ss in 1:ngibbs){
    lp[ss] <- likelihood(gibbsSample$z_gibbs[,,thinInd[ss]], gibbsSample$R_gibbs[,,thinInd[ss] ] )$log.lik
    if(ss%%2500==0){print(ss)}
  }
  lp5[ch,] <- lp
  pi_gibbs[,,(ch-1)*ngibbs+ 1:ngibbs] <- gibbsSample$E_gibbs[,,thinInd]
  sig2_gibbs[(ch-1)*ngibbs+ 1:ngibbs] <- gibbsSample$sig2_gibbs[thinInd]
  C_gibbs[,,(ch-1)*ngibbs+ 1:ngibbs] <- gibbsSample$C_gibbs[,,thinInd]
  
  T_post_mean[[ch]] <- apply( gibbsSample$U_gibbs[,,thinInd],c(1,2),mean)
  rm(gibbsSample) 
}


# Color for traceplotß
tracecol <- c( rgb(0.1,0.1,0.1,0.7), rgb(0.8,0.1,0.1,0.5),
               rgb(0.1,0.8,0.1,0.7), rgb(0.1,0.1,0.8,0.5))

#pdf("tracePlot4.pdf",width=6,height=4)
plot( lp5[1,],type="l",lty=1, col=tracecol[1],
      main="Trace plot of", xlab="", ylab="", ylim=c(-6750,-6300))
lines( lp5[2,],type="l",lty=1,col=tracecol[2])
lines( lp5[3,],type="l",lty=1,col=tracecol[3])
lines( lp5[4,],type="l",lty=1,col=tracecol[4])
legend("topright",paste("Chain",1:4), lty=1,col=tracecol, lwd=2)
#dev.off()




## Posterior means of edge inclusion probabilities (pi_jk)
pi_hat   <- Reduce("+", alply(pi_gibbs,3) )/(ngibbs*nchain)
utriFlag <- upper.tri(pi_hat) # Flag for upper triangle part components
ltriFlag <- lower.tri(pi_hat) # Flag for lower triangle part components
pi_vec   <- pi_hat[utriFlag]

## Roughly bimodal concentrated near 0 and 1.
hist( pi_vec,50)



## False discovery rate (Mitra et al., 2013; Peterson et al., 2015)
fdr_th <- function(pi_vec,th){
	num   <- sum( (1-pi_vec)*( pi_vec>th) )
	denom <- sum( pi_vec>th ) #+ 1e-8
	fdr <- num/denom
	return(fdr)
}


ths <- seq(min(pi_vec),max(pi_vec),l=500)
plot(ths, Vectorize(function(x) fdr_th(pi_vec,x)) (ths), type="l" )
abline(h=0.1,col=2)


alpha_target <- 0.1 #Target FDR
indexFDRlessthan0.1 <- which( Vectorize(function(x) fdr_th(pi_vec,x)) (ths) <= alpha_target)
thresh <- min(ths[indexFDRlessthan0.1]) # The minimum graph threshold with FDR<= alpha_target
thresh
round(thresh,3)
## Check
## Shoud be negative to control FDR below alpha_target
fdr_th(pi_vec,thresh) - alpha_target


## Estimate of adjacency matrix 
Wtmp  <- 1*( (Reduce("+" , alply(pi_gibbs,3) )/(ngibbs*nchain) )> thresh )
mean(Wtmp[utriFlag])
What <- Wtmp
diag(What) <- 0
colnames(What) <- colnames(x)
rownames(What) <- colnames(x)  
grph <- graph_from_adjacency_matrix(What, mode="undirected")


# Try various graph layouts of the estimated network
set.seed(77845)
ceb <- cluster_edge_betweenness(grph) 

## Draw graph with layout = layout_with_kk
V(grph)$label.cex <- 1.1
V(grph)$label.color <- "black"


#which(colnames(What)=="Methanosphaera")
#V(grph)$label.dist[54] <- -1

#Change the label location so that all fit genera in 
V(grph)$label.dist <- rep(1.2,54)
V(grph)$edge.color <- "grey"

layout_kk <- layout_with_kk(grph)


layout_kk[,2] <- -layout_kk[,2]

## Community memberships in phylogenetic tree
cols <- membership(ceb)
cols[ match( ceb[[1]] , qmptree$tip.label ) ] <- "black"
cols[ match( ceb[[2]] , qmptree$tip.label ) ] <- "cyan2"
cols[ match( ceb[[3]] , qmptree$tip.label ) ] <- "black"
cols[ match( ceb[[4]] , qmptree$tip.label ) ] <- "chartreuse2"
cols[ match( ceb[[5]] , qmptree$tip.label ) ] <- "red"
cols[ match( ceb[[6]] , qmptree$tip.label ) ] <- "black"

shapes <- rep(1,p)
shapes[ match( ceb[[1]] , qmptree$tip.label ) ] <- 3
shapes[ match( ceb[[3]] , qmptree$tip.label ) ] <- 3
shapes[ match( ceb[[6]] , qmptree$tip.label ) ] <- 3


V(grph)$shape <- vertex.shapes()[shapes]
#save(layout_kk,file="graphLayout.RData")
#pdf(paste(figPath,"network.tree.4chains.pdf",sep=""),width=12.5,height=5)
par(oma=c(1,1,1,1),mar=c(0,0,1,1))
plot.igraph(grph, vertex.label=colnames(What), vetex.label.color="black", 
     vertex.color = cols,
     asp=0.45, edge.width=0.3, layout = layout_kk , vertex.size=rep(3,54), main="PhyloBCG", rescale=TRUE) 
dev.off()


# This is for a presentation slide            
#pdf(paste(figPath,"nomembership.tree.4chains.pdf",sep=""),width=14,height=9)
par(mar=c(1,1,1,0),xpd=TRUE,oma=c(1,2,1,0))
  h <- plot(qmptree, type="phylogram",tip.col="black",x.lim=9.5)
  ape::nodelabels(qmptree$node.label, cex=0.8, adj=1)
#dev.off()


## Community memberships in phylogenetic tree
  cols <- rep("black",length(qmptree$tip.label))
  cols[ match( ceb[[2]] , qmptree$tip.label ) ] <- "cyan2"
  cols[ match( ceb[[4]] , qmptree$tip.label ) ] <- "chartreuse2"
  cols[ match( ceb[[5]] , qmptree$tip.label ) ] <- "red"
  

  mem <- rep("*",length(qmptree$tip.label))
  mem[ match( ceb[[2]] , qmptree$tip.label ) ] <- "(2)"
  mem[ match( ceb[[4]] , qmptree$tip.label ) ] <- "(1)"
  mem[ match( ceb[[5]] , qmptree$tip.label ) ] <- "(3)"

  
# Label community membership  
labeled.qmptree <- qmptree
labeled.qmptree$tip.label <- paste(mem, qmptree$tip.label)
labeled.qmptree$tip.label <- paste(mem, qmptree$tip.label)
#pdf(paste(figPath,"membership.tree.4chains.pdf",sep=""),width=16,height=12)
par(mar=c(1,1,1,0),xpd=TRUE,oma=c(1,2,1,0))
  h <- plot(labeled.qmptree, type="phylogram",tip.col=cols,x.lim=9.5)
  ape::nodelabels(labeled.qmptree$node.label, cex=0.8, adj=1)
#dev.off()






###########################################################################
############### Draw posterior latent positions
###########################################################################

shapes <- rep(15,length(qmptree$tip.label))
shapes[ match( ceb[[2]] , qmptree$tip.label ) ] <- 16
shapes[ match( ceb[[4]] , qmptree$tip.label ) ] <- 17
shapes[ match( ceb[[5]] , qmptree$tip.label ) ] <- 18


post_latentpositions <- Reduce("+", T_post_mean )/nchain
colnames( post_latentpositions ) <- colnames(x)
tmp <- data.frame( x1 = post_latentpositions[1,], x2 = post_latentpositions[2,] )


library(ggrepel)
library(ggfortify)#install.packages('ggfortify')

nudge_x <-  nudge_y <- rep(0,p)

nudge_x[cols == "red"] <-  -1.2
nudge_y[cols == "red"] <-  -1.2

nudge_x[cols == "cyan2"] <-   2.0
nudge_y[cols == "cyan2"] <-  -1.0

nudge_y[cols == "chartreuse2"] <-   2.2

post.position <-  ggplot(tmp,aes(x=x1, y=x2, fill=cols))# + stat_ellipse(type = "norm", linetype = 2) 
latent.positions <- post.position +
                    geom_point(aes(x=x1,y=x2), color=cols, shape=shapes, size=3) +
                                ylim(lower=-4,upper=5.0) + xlim(lower=-5,upper=5.5) +
                    geom_label_repel(aes(x=x1,y=x2,label=colnames(x) ), size=3, 
                                    nudge_y = nudge_y,  nudge_x = nudge_x, 
                                    arrow = arrow(length = unit(0.01, "npc")),
                                    box.padding = 0.01, point.padding = 0.2, label.size=0.2, 
                                    segment.color=rep("grey50",p), 
                                    max.time=50, alpha=0.5, max.overlaps = 100) +
                   xlab("") + ylab("") +
                   ggtitle("Posterior mean latent positions")  +
                   theme( plot.title = element_text(hjust = 0.5) )
latent.positions 

#ggsave(paste(figPath,"latent_positions.pdf",sep=""),latent.positions, width=10,height=6)




###########################################################################
############### Draw heatmap of posterior partial correlation matrix
###########################################################################
OmegaHat <-  Reduce("+" , alply(C_gibbs,3) )/(ngibbs*nchain)


# Keep the diagonal elements
diagtmp <- diag(OmegaHat)
# Make the precision matrix sparse by multiplying the adjacency matrix What
OmegaHat <- OmegaHat*What
# diag(What) = zeros
# restore diagonal elements
diag(OmegaHat) <- diagtmp

# Again set diagonal elements of the partial correlation matrix 0s
parCor <- -cov2cor(OmegaHat)
diag(parCor) <- 0

rownames(parCor) <- colnames(x) # column and row names
colnames(parCor) <- colnames(x)
save(parCor, file = paste(resultsPath,"pcor.tree_4chains.RData",sep=""))








# create a coordinate vectors for geom_tile
TileMat <- parCor
coord <- cbind(rowv = rep(rownames(parCor), each = p), 
               colv = rep(rownames(parCor), times = p))
# Save vectorized partial correlation matrix to data frame with corresponding coordinates
df_Tilemat <- cbind.data.frame(coord, parC = c(parCor))

# rows start from bottom to top and cols start from left to right.
# but we want start from the left top corner for both triangles. So make colv factor in reverse level.
df_Tilemat$rowv <- factor(df_Tilemat$rowv, levels = rownames(parCor))
df_Tilemat$colv <- factor(df_Tilemat$colv, levels = rev(rownames(parCor)))


## global aes
tile.pc.hat <- ggplot(df_Tilemat, aes(x = rowv, y = colv), fill = parCor)
## to get the rect filled
tile.pc.hat <- tile.pc.hat + geom_tile(aes(fill = parCor), color = "black", alpha = 1, width = 1, height = 1)
tile.pc.hat <- tile.pc.hat + scale_fill_gradient2(high = "blue", low = "red", name = "ParCor", limits=c(-.5,.4)) 
## color of the corresponding aes
tile.pc.hat <- tile.pc.hat +   coord_cartesian(xlim = c(1, p), ylim = c(1, p), clip = 'off') 
tile.pc.hat <- tile.pc.hat +   theme_bw() ## Add boundary
## Add legend
tile.pc.hat <- tile.pc.hat + theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0), 
                                   axis.title.x=element_blank(),
                                   axis.title.y=element_blank(),
                                   legend.position = c(0, 1.2), 
                                   legend.justification="left", 
                                   legend.direction="horizontal", 
                                   plot.margin = margin(1, 1, 1, 0, "cm")) 
#annotate("text", x = (p+2)/2, y = -0.7, label = "BCG on Quantitative counts") + 
tile.pc.hat <- tile.pc.hat + annotate("text", x = (p+2)/2, y = -0.7, angle = 0, label = "Partial Correlation")
tile.pc.hat <- tile.pc.hat + scale_x_discrete(position = "top") # Move x tick labels to the top
tile.pc.hat <- tile.pc.hat + scale_alpha_continuous(guide=FALSE) +  scale_size_continuous(guide=FALSE) #Better graphical representation
tile.pc.hat <- tile.pc.hat + geom_segment(aes(x=0.5, y=(p+0.5), xend = (p+0.5), yend=0.5), color="grey", size = 0.5) #Draw diagonal line
tile.pc.hat

#ggsave(paste(figPath,"parCor.tree_1e-3.pdf",sep=""),tile.pc.hat, width=12,height=9)



###########################################################################
############### Draw heatmap of posterior edge inclusion probability
###########################################################################
pi_hat <- Reduce("+", alply(gibbsSample$pi_gibbs[,,],3) )/nmc
diag(pi_hat) <- 0

rownames(pi_hat) <- colnames(x) # column and row names
colnames(pi_hat) <- colnames(x)


# create a coordinate vectors for geom_tile
TileMat <- pi_hat
coord <- cbind(rowv = rep(rownames(pi_hat), each = p), 
               colv = rep(rownames(pi_hat), times = p))
# Save vectorized pi_hat matrix to data frame with corresponding coordinates
df_Tilemat <- cbind.data.frame(coord, pi_hat = c(pi_hat))

# rows start from bottom to top and cols start from left to right.
# but we want start from the left top corner for both triangles. So make colv factor in reverse level.
df_Tilemat$rowv <- factor(df_Tilemat$rowv, levels = rownames(pi_hat))
df_Tilemat$colv <- factor(df_Tilemat$colv, levels = rev(rownames(pi_hat)))



tile.pi.hat <- ggplot(df_Tilemat, aes(x = rowv, y = colv), fill = pi_hat)         ## global aes
tile.pi.hat <- tile.pi.hat + geom_tile(aes(fill = pi_hat), color = "black", alpha = 1, width = 1, height = 1)## to get the rect filled
tile.pi.hat <- tile.pi.hat + scale_fill_gradient2(high = "blue", low = "white", name = expression(pi[jk]), limits=c(0,1)) 
  ## color of the corresponding aes
tile.pi.hat <- tile.pi.hat +   coord_cartesian(xlim = c(1, p), ylim = c(1, p), clip = 'off') 
tile.pi.hat <- tile.pi.hat +   theme_bw() ## Add boundary
## Add legend
tile.pi.hat <- tile.pi.hat + theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = c(0, 1.2), 
        legend.justification="left", 
        legend.direction="horizontal", 
        plot.margin = margin(1, 1, 1, 0, "cm")) 
  #annotate("text", x = (p+2)/2, y = -0.7, label = "BCG on Quantitative counts") + 
tile.pi.hat <- tile.pi.hat + annotate("text", x = (p+2)/2, y = -0.7, angle = 0, label = "Posterior edge inclusion probability")
tile.pi.hat <- tile.pi.hat + scale_x_discrete(position = "top")
tile.pi.hat <- tile.pi.hat + scale_alpha_continuous(guide=FALSE) +  scale_size_continuous(guide=FALSE)
tile.pi.hat <- tile.pi.hat + geom_segment(aes(x=0.5, y=(p+0.5), xend = (p+0.5), yend=0.5), color="grey", size = 0.5)
tile.pi.hat

#ggsave(paste(figPath,"post.pi_jk.pdf",sep=""),tile.pi.hat, width=12,height=9)


#ggsave(paste(figPath,"tree.pi.parCor.pdf",sep=""), grid.arrange(tile.pi.hat,tile.pc.hat, ncol=2), width=18,height=9)









###########################################################################
############### Draw heatmap of posterior covariance matrix
###########################################################################
Sighat <- Reduce("+", alply(gibbsSample$R_gibbs[,,],3) )/nmc

rownames(Sighat) <- colnames(x) # column and row names
colnames(Sighat) <- colnames(x)


# create a coordinate vectors for geom_tile
TileMat <- Sighat
coord <- cbind(rowv = rep(rownames(Sighat), each = p), 
               colv = rep(rownames(Sighat), times = p))
# Save vectorized pi_hat matrix to data frame with corresponding coordinates
df_Tilemat <- cbind.data.frame(coord, Sighat = c(Sighat))

# rows start from bottom to top and cols start from left to right.
# but we want start from the left top corner for both triangles. So make colv factor in reverse level.
df_Tilemat$rowv <- factor(df_Tilemat$rowv, levels = rownames(Sighat))
df_Tilemat$colv <- factor(df_Tilemat$colv, levels = rev(rownames(Sighat)))



tile.S.hat <- ggplot(df_Tilemat, aes(x = rowv, y = colv), fill = Sighat)         ## global aes
tile.S.hat <- tile.S.hat + geom_tile(aes(fill = Sighat), color = "black", alpha = 1, width = 1, height = 1)## to get the rect filled
tile.S.hat <- tile.S.hat + scale_fill_gradient2(high = "blue", low = "red", name = expression(sigma[jk]), limits=c(-1,1)) 
## color of the corresponding aes
tile.S.hat <- tile.S.hat +   coord_cartesian(xlim = c(1, p), ylim = c(1, p), clip = 'off') 
tile.S.hat <- tile.S.hat +   theme_bw() ## Add boundary
## Add legend
tile.S.hat <- tile.S.hat + theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0), 
                                   axis.title.x=element_blank(),
                                   axis.title.y=element_blank(),
                                   legend.position = c(0, 1.2), 
                                   legend.justification="left", 
                                   legend.direction="horizontal", 
                                   plot.margin = margin(1, 1, 1, 0, "cm")) 
#annotate("text", x = (p+2)/2, y = -0.7, label = "BCG on Quantitative counts") + 
tile.S.hat <- tile.S.hat + annotate("text", x = (p+2)/2, y = -0.7, angle = 0, label = "Posterior correlation")
tile.S.hat <- tile.S.hat + scale_x_discrete(position = "top")
tile.S.hat <- tile.S.hat + scale_alpha_continuous(guide=FALSE) +  scale_size_continuous(guide=FALSE)
tile.S.hat <- tile.S.hat + geom_segment(aes(x=0.5, y=(p+0.5), xend = (p+0.5), yend=0.5), color="grey", size = 0.5)
tile.S.hat





