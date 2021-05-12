rm(list=ls())
#setwd("/Users/heecheolchung/Dropbox/Research/TAMU/graphical/p50cdf/")
#funcPath   <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/functions/"
#resultPath <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/p50cdf/"

setwd("/general/home/hcchung/graphical/p50cdf/m5/")
funcPath   <- "/general/home/hcchung/graphical/functions/"
resultPath <- "/general/home/hcchung/graphical/p50cdf/results/"

##Load required libraries
source(paste(funcPath,"libraries.R",sep=""))




source(paste(funcPath,"genTree.R",sep=""))
nlam <- 100
lambda.max <- 0.9
lambda.min <- 0.001
lambdaseq <- exp(seq(log(lambda.max), log(lambda.min), length = nlam))


nrep <- 50 # Number of replicated datasets


set.seed(77843,kind = "Mersenne-Twister" ,sample.kind = "Rejection" )
ndata <- 10
seeds <- sample(19999:29999,ndata)#print(seeds)# [1] 23804 25990 23816 27180 25273 26819 20145 26898 25655 22684
print(seeds)

cscale <- 3 # Tree scale parameter

## Load QMP data for synthetic data generation
# dim(QMP) #There are 106 subjects and 54 genera
load("QMPtree.RData")
## Variables with to small zero proportions -> make they have at least 20%
## Reason : If one has too small zero proportion then the threshold value will be also very small
##          This small threshold value causes numerical problem in sampling from mv truncated normal.
##          Ex. N(0,1) truncated above at -4, i.e., x \in (-Inf, -4)
# This is the source code that reduced QMP data dimension from 54 to 50
source(paste(funcPath,"QMPreduced_for_sim.R",sep=""))

n <- dim(QMP)[1] # Sample size
p <- dim(QMP)[2]  # Data dimension reduced from 54 to 50
K <- 2   # Latent space dimension



iseed=1;
print(c(iseed,seeds[iseed]))
set.seed( seeds[iseed], kind = "Mersenne-Twister" ,sample.kind = "Rejection" )
source(paste(funcPath,"TreeDataGeneration.R",sep=""))
#plot(graph_from_adjacency_matrix( Wtrue, mode="undirected"))


bt <- Sys.time()
########################################################################
########           Data Generation                  ####################
########################################################################

x50 <- vector("list",nrep)
for( ii in 1:nrep){
  x50[[ii]] <- synthData_from_ecdf(QMP, mar = 2, SigmaTrue , n=n, seed = NULL, verbose = FALSE)
}
# Fit SPRING for each 50 replicated data sets
Rep50 <- foreach( irep = 1:nrep)%dopar%{

        x <- x50[[irep]]
        fit.spring <- SPRING(x, Rmethod = "approx", quantitative = TRUE, verbose = FALSE,
                            lambdaseq = "data-specific", nlambda = 100 )

        # StARS-selected lambda index based on the threshold (default = 0.1)
        opt.K <- fit.spring$output$stars$opt.index
        # Estimated adjacency matrix from sparse graphical modeling technique ("mb" method) (1 = edge, 0 = no edge)
        adj.K <- as.matrix(fit.spring$fit$est$path[[opt.K]])
        # Estimated adjacency matrix
        Wtmp <- adj.K
        # Upper triangluar part of the estimated graph
        edgeHat  <- Wtmp[ upper.tri(Wtmp) ]
        # Upper triangluar part of the true graph
        edgeTrue <- as.matrix(Wtrue)[ upper.tri(Wtmp) ]
        # The total number of connected edges
        edgeTotal <- sum( edgeTrue==1 )

        TP  <- sum( edgeHat[ edgeTrue==1 ] ) # of true positives
        FP  <- sum( edgeHat[ edgeTrue==0 ] ) # of false positives
        MCC <- cor(edgeHat,edgeTrue) # Mattew correlation coefficient
        

        What <- Wtmp

        simSummary <- list(
                      fit.spring = fit.spring,
                      What = What,
                      Wtrue = Wtrue,
                      edgeHat = edgeHat,
                      edgeTrue = edgeTrue,
                      edgeTotal = edgeTotal,
                      TP = TP,
                      FP = FP,
                      MCC = MCC
                        )

            return( simSummary )

        }
	et <- Sys.time()

print(et-bt)

save.image(paste(resultPath,"SPRING_",iseed,".RData",sep=""))






