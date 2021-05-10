rm(list=ls())
library(abind)#install.packages("abind")
library(plyr)
library(xtable)
library(gridExtra)
setwd("/Users/heecheolchung/Dropbox/Research/TAMU/graphical/p50cdf/results/")


## False discovery rate (Mitra et al., 2013; Peterson et al., 2015)
fdr_th <- function(pi_vec,th){
  num   <- sum( (1-pi_vec)*( pi_vec>th) )
  denom <- sum( pi_vec>th ) + 1e-8
  fdr <- num/denom
  return(fdr)
}



##################################
##     Simulation Summary       ##
##################################
allnorms <- function(d) vapply(c("1","I","F","M","2"),
                               norm, x = d, double(1))

mods <- c("O_","D_","T_")# Oracle, Dist, Tree
settingNames <- c( paste("Gibbs_cdf",mods,sep=""),  "SPRING_", "SPIECEASI_")
#settingNames <- c( paste("Gibbs_cdf",mods,sep=""), "SPRING_" )

ndata        <- 10

MCC   <- matrix(0,ndata,length(settingNames))# Save 10 average MCC values
MCCse <- matrix(0,ndata,length(settingNames))

TPR   <- matrix(0,ndata,length(settingNames))
TPRse <- matrix(0,ndata,length(settingNames))

FPR   <- matrix(0,ndata,length(settingNames))
FPRse <- matrix(0,ndata,length(settingNames))

FDR   <- matrix(0,ndata,length(settingNames))
FDRse <- matrix(0,ndata,length(settingNames))

p <- 50
nUpperTri <- p*(p-1)/2
upperFlag <- upper.tri(matrix(NA,p,p))
allsumm <- matrix(0,ndata, length(settingNames) + 1)

colnames(allsumm) <- c( rep(c("Fixed","Dist","Tree", "SPRING", "SPIECEASI"), each=1 ), "Marg.Inc.Prob" )

mean.sig2 <-  rep(0,10)
mean.gamma <- rep(0,10)

for( iset in 1:length(settingNames) ){

	for( idata in 1:ndata ){


		load( paste(settingNames[iset],idata,".RData",sep="") )

		if( iset == 4 | iset == 5 ){ # SPEIC-EASI (4) or SPRING (5)
 		  MCC[idata,iset]   <- mean( unlist( lapply( Rep50, function(x) x$MCC ) ), na.rm=TRUE )
 		  MCCse[idata,iset] <- sd( unlist( lapply( Rep50, function(x) x$MCC ) ), na.rm=TRUE  )/sqrt(nrep)
 		
 		  TPR[idata,iset]   <- mean( unlist( lapply( Rep50, function(x) x$TP/x$edgeTotal  ) ) )
 		  TPRse[idata,iset] <- sd( unlist( lapply( Rep50, function(x) x$TP/x$edgeTotal  ) ) )/sqrt(nrep)
 		
 		  FPR[idata,iset]   <- mean( unlist( lapply( Rep50, function(x) x$FP/(nUpperTri - x$edgeTotal)  ) ) )
 		  FPRse[idata,iset] <- sd( unlist( lapply( Rep50, function(x) x$FP/(nUpperTri - x$edgeTotal)  ) ) )/sqrt(nrep)

 		  FDR[idata,iset]   <- mean( unlist( lapply( Rep50, function(x) x$FP/(x$TP + x$FP)  ) ) )
 		  FDRse[idata,iset] <- sd( unlist( lapply( Rep50, function(x) x$FP/(x$TP + x$FP)  ) ) )/sqrt(nrep)
 		  
		}else{
    #######################################################################################
		  alpha_target <- 0.05
		  
		  pihatList <- lapply(Rep50, function(x) x$What[upperFlag] )
		  ths  <- seq(0.01,0.99,l=100)
      fdrc <- rep(0,50)
      for( pp in 1:nrep ){
      indexFDRlessthan0.05     <- which( Vectorize(function(x) fdr_th(pihatList[[pp]],x)) (ths) <= alpha_target)
      fdrc[[pp]]               <- min(ths[indexFDRlessthan0.05])
      }

    graphs <- mapply(function(x,y) x>y, lapply(Rep50, function(x) x$What ),alply(fdrc,1), SIMPLIFY = FALSE)
    #######################################################################################
    
	  trueW  <- lapply(Rep50, function(x) x$Wtrue )
	  	  
	  mccs <- mapply( function(x,y) cor(x[upperFlag], as.matrix(y)[upperFlag]), graphs,  trueW)
	  tprs <- mapply( function(x,y) sum(x[upperFlag][ as.matrix(y)[upperFlag]!=0])/sum(as.matrix(y)[upperFlag]!=0), 
	                  graphs,  trueW)
	  fprs <- mapply( function(x,y) sum(x[upperFlag][ as.matrix(y)[upperFlag]==0])/sum(as.matrix(y)[upperFlag]==0), 
	                  graphs,  trueW)

	  discv.t <- mapply( function(x,y) sum(x[upperFlag][ as.matrix(y)[upperFlag]!=0]), graphs,  trueW)
	  discv.f <- mapply( function(x,y) sum(x[upperFlag][ as.matrix(y)[upperFlag]==0]), graphs,  trueW)
	  	  
	  MCC[idata,iset]   <- mean( mccs )
	  MCCse[idata,iset] <- sd( mccs  )/sqrt(nrep)
	  
	  TPR[idata,iset]   <- mean( tprs )
	  TPRse[idata,iset] <- sd( tprs )/sqrt(nrep)
	  
	  FPR[idata,iset]   <- mean( fprs )
	  FPRse[idata,iset] <- sd( fprs )/sqrt(nrep)

	  FDR[idata,iset]   <- mean( discv.f/(discv.t+discv.f) )
	  FDRse[idata,iset] <- sd( discv.f/(discv.t+discv.f) )/sqrt(nrep)
	  
	  if(iset == 3 ){
	    mean.sig2[idata] <- mean( unlist( lapply(Rep50, function(x) mean(x$sig2_gibbs)) ) )
	  }else if(iset ==2 ){
	    mean.gamma[idata] <- mean( unlist( lapply(Rep50, function(x) x$gamma[3] ) ) )
	    }
	  
		}

		tmp <- format( round( c(MCC[idata,iset], TPR[idata,iset], FPR[idata,iset]), 3), nsmall=3 )
		tmp <- paste(tmp[1]," (",tmp[2],"/",tmp[3],")",sep="")
		
		allsumm[idata,iset] <- tmp
		
		if( iset == 1){# iset can be any setting index
		allsumm[idata,length(settingNames)+1]    <- format( round( Rep50[[1]]$edgeTotal/nUpperTri, 3), nsmall=3)
		}

		}
		
}	

## Extract latex table
xtable( allsumm, digits=3 )

## Load clustering coefficients for reodering
CC <- read.csv("/Users/heecheolchung/Dropbox/Research/TAMU/graphical/p50cdf/clusteringCoefficient.csv")[,1]


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols <- gg_color_hue(9)
cols <- c(cols[1],cols[4],cols[7])

meths    <- c("Oracle","Dist","PhyloBCG","SPRING","SPIEC-EASI")
tree.order   <- order(CC,decreasing=TRUE)




## Draw summary plot
library(ggplot2)
library(tidyr)


MCC   <- MCC[tree.order,]
MCCse <- MCCse[tree.order,]

TPR   <- TPR[tree.order,]
TPRse <- TPRse[tree.order,]

FPR   <- FPR[tree.order,]
FPRse <- FPRse[tree.order,]

FDR   <- FDR[tree.order,]
FDRse <- FDRse[tree.order,]


df.1 <- data.frame( mcc =as.vector(MCC), mccse = 2*as.vector(MCCse), 
                    tpr =as.vector(TPR), tprse = 2*as.vector(TPRse), 
                    fpr =as.vector(FPR), fprse = 2*as.vector(FPRse),
                    fdr =as.vector(FDR), fdrse = 2*as.vector(FDRse),
                    Model=rep(meths,each=10), tree.id=rep(1:10,time=length(meths)) )
df.1$Model <- factor(df.1$Model, level = meths[c(3,2,1,4,5)])


# dev.off()
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#Summary of average MCC
ggmcc   <- ggplot(df.1, aes(x=tree.id, y=mcc, group=Model, color=Model)) + 
  geom_line()  +
  geom_point(size=0.5) +
  geom_errorbar(aes(ymin=mcc-mccse, ymax=mcc+mccse), width=.5,
                position=position_dodge(0.01))

ggmcc <- ggmcc + 
  labs(x = "Phylogenetic tree index", y = "MCC") + 
  ylim(low=0.01,high=0.99) +
  #xlim(low=0.00,high=12) +
  scale_x_continuous( breaks=1:10) +
  ggtitle("Average Matthews correlation coefficient") +
  theme(plot.title = element_text(hjust=0.5), legend.position = c(0.9,0.8))+ 
  guides(color = guide_legend(override.aes = list(size = 1.2),title=""))  # Legned line thickness
# if limits = c(1,10), then the 10th stand error bars are truncated so I give little more space 
# by setting c(1,10.2)
#----------------------------------------------------------------------------------------
#Summary of average true positive rates
ggtpr   <- ggplot(df.1, aes(x=tree.id, y=tpr, group=Model, color=Model)) + 
  geom_line()  +
  geom_point(size=0.5) +
  geom_errorbar(aes(ymin=tpr-tprse, ymax=tpr+tprse), width=.5,
                position=position_dodge(0.01))

ggtpr <- ggtpr + 
  labs(x = "Phylogenetic tree index", y = "TPR") + 
  ylim(low=0.00,high=1) +
  scale_x_continuous( breaks=1:10) +
  ggtitle("Average true positive rate") +
  theme(plot.title = element_text(hjust=0.5), legend.position = c(1.6,0.8)) 

#----------------------------------------------------------------------------------------
#Summary of average false positive rates
ggfpr   <- ggplot(df.1, aes(x=tree.id, y=fpr, group=Model, color=Model)) + 
  geom_line()  +
  geom_point(size=0.2) +
  geom_errorbar(aes(ymin=fpr-fprse, ymax=fpr+fprse), width=.5,
                position=position_dodge(0.01))

ggfpr <- ggfpr + 
  labs(x = "Phylogenetic tree index", y = "FPR") + 
  ylim(low=0.00,high=0.2) +
  scale_x_continuous( breaks=1:10) +
  ggtitle("Average false positive rate") +
  theme(plot.title = element_text(hjust=0.5), legend.position = c(1.6,0.5))

#----------------------------------------------------------------------------------------
#Summary of average false discovery rates
ggfdr   <- ggplot(df.1, aes(x=tree.id, y=fdr, group=Model, color=Model)) + 
  geom_line()  +
  geom_point(size=0.2) +
  geom_errorbar(aes(ymin=fdr-fdrse, ymax=fdr+fdrse), width=.5,
                position=position_dodge(0.01))

ggfdr <- ggfdr + 
  labs(x = "Phylogenetic tree index", y = "FDR") + 
  ylim(low=0.00,high=0.5) +
  scale_x_continuous( breaks=1:10) +
  ggtitle("Average false discovery rate") +
  theme(plot.title = element_text(hjust=0.5), legend.position = c(1.6,0.5)) + 
  geom_abline(intercept=0.05,slope=0, linetype = "dashed")
#----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

#grid.arrange(ggmcc,ggtpr,ggfdr,ggfpr, ncol=2)
#grid.arrange(arrangeGrob(ggmcc,nrow=1,ncol=1),arrangeGrob(ggtpr,ggfpr,nrow=2,ncol=1), arrangeGrob(ggfdr,nrow=1,ncol=1), widths=c(1.5,1) )
#grid.arrange(arrangeGrob(ggmcc,nrow=1,ncol=1),arrangeGrob(ggtpr,ggfpr,nrow=2,ncol=1), widths=c(1.5,1) )
#grid.arrange(arrangeGrob(ggmcc,nrow=2,ncol=1),arrangeGrob(ggtpr,ggfpr,nrow=1,ncol=2), widths=c(1.5,1) ,ncol=2)
lay <- rbind(c(1,1),
             c(2,3) )
grid.arrange(ggmcc,ggtpr,ggfdr, layout_matrix=lay, ncol=2)



#pdf("simsum_cdf_all_ejk_fdr_5.pdf",height = 7, width=12)
lay <- rbind(c(1,1),
             c(2,3) )
grid.arrange(ggmcc,ggtpr,ggfpr, layout_matrix=lay, ncol=2)
  #dev.off()




getwd()











#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# For presentation slides
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
slidetree <- c(2,9)

df.1 <- data.frame( mcc =as.vector(MCC[slidetree,1:3]), mccse = 2*as.vector(MCCse[slidetree,1:3]), 
                    tpr =as.vector(TPR[slidetree,1:3]), tprse = 2*as.vector(TPRse[slidetree,1:3]), 
                    fpr =as.vector(FPR[slidetree,1:3]), fprse = 2*as.vector(FPRse[slidetree,1:3]),
                    fdr =as.vector(FDR[slidetree,1:3]), fdrse = 2*as.vector(FDRse[slidetree,1:3]),
                    Model=rep(meths[1:3],each=2), tree.id=rep(1:2,time=length(meths[1:3])) )
df.1$Model <- factor(df.1$Model, level = meths[c(3,2,1)])


# dev.off()
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#Summary of average MCC
ggmcc   <- ggplot(df.1, aes(x=tree.id, y=mcc, group=Model, color=Model)) + 
  #geom_line()  +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mcc-mccse, ymax=mcc+mccse), width=.2)

ggmcc <- ggmcc + 
  labs(x = "Phylogenetic tree index", y = "MCC") + 
  ylim(low=0.01,high=0.99) +
  #xlim(low=0.00,high=12) +
  scale_x_continuous( breaks=1:10) +
  #ggtitle("Average Matthews correlation coefficient") +
  theme(plot.title = element_text(hjust=0.5), legend.position = c(0.8,0.9))+ 
  guides(color = guide_legend(override.aes = list(size = 1.2),title=""))  # Legned line thickness
# if limits = c(1,10), then the 10th stand error bars are truncated so I give little more space 
# by setting c(1,10.2)
#----------------------------------------------------------------------------------------
#Summary of average true positive rates
ggtpr   <- ggplot(df.1, aes(x=tree.id, y=tpr, group=Model, color=Model)) + 
  geom_point(size=2) +
  geom_errorbar(aes(ymin=tpr-tprse, ymax=tpr+tprse), width=.2)

ggtpr <- ggtpr + 
  labs(x = "Phylogenetic tree index", y = "TPR") + 
  ylim(low=0.00,high=1) +
  scale_x_continuous( breaks=1:10) +
  #ggtitle("Average true positive rate") +
  theme(plot.title = element_text(hjust=0.5), legend.position = c(1.6,0.8)) 

#----------------------------------------------------------------------------------------
#Summary of average false positive rates
ggfpr   <- ggplot(df.1, aes(x=tree.id, y=fpr, group=Model, color=Model)) + 
  geom_point(size=2) +
  geom_errorbar(aes(ymin=fpr-fprse, ymax=fpr+fprse), width=.2)

ggfpr <- ggfpr + 
  labs(x = "Phylogenetic tree index", y = "FPR") + 
  ylim(low=0.00,high=1) +
  scale_x_continuous( breaks=1:10) +
  #ggtitle("Average false positive rate") +
  theme(plot.title = element_text(hjust=0.5), legend.position = c(1.6,0.5))



#pdf("simsum_cdf_tree_2_9.pdf",height = 6, width=10)
grid.arrange(ggmcc,ggtpr,ggfpr, ncol=3)
#dev.off()

getwd()
