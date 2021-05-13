rm(list=ls())
#setwd("/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/")
#funcPath   <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/functions/"
#resultPath <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/QMPdata/results/"

setwd("//general/home/hcchung/graphical/adaQMP/")
funcPath   <- "/general/home/hcchung/graphical/functions/"
resultPath <- "/general/home/hcchung/graphical/adaQMP/results/"

library(boot); library(tmvtnorm); #install.packages("tmvtnorm")
library(huge) #install.packages("huge")
library(igraph); library(BDgraph) #install.packages("BDgraph") #install.packages("igraph")
library(plyr);library(ggplot2);library(reshape);library(gridExtra);
library(truncnorm)
# Load functions
source(paste(funcPath,"conditionalSample.R",sep=""))
source(paste(funcPath,"SSVS.R",sep=""))
source(paste(funcPath,"AlbertChib_mrg.R",sep=""))
source(paste(funcPath,"treeGibbs.R",sep=""))
# Load QMP data 
load("QMPtree.RData")

x  <- QMP
np <- dim(x)
n  <- np[1] # Sample size
p  <- np[2]  # Data dimension
K <- 2   # Latent space dimension


H_t     <- cov2cor(H) # Tree correlation matrix
H_tinv  <- solve(H_t) # Tree precision matrix


########################################################################
########         Hyper parameters (fixed)               ################
########################################################################
# Fixed hyperparameters
# Note that the parameterization used in the code is slightly different from those in Wang (2014).  )
h = 2500    # (h in code) = (h in paper)^2
lambda = 1    # Prior for diagonal entries of the precision matrix wjj ~ Exp(wjj|lambda/2)
IGsig2 = rep(1e-3,2) # Inv-gamma parameters for tree scale parameter (sigma2)
IGv0   = rep(1e-3,2) # Inv-gamma parameters for the spike and slab (v0)
hyperparameters = list(h=h, lambda=lambda, IGsig2=IGsig2, IGv0=IGv0)


########################################################################
########    Pseudo Data Generation (zhat)          ####################
########################################################################

delta_mc <- qnorm( colSums( x == 0 )/n )

# Get the scaled, ( n/(n+1) ), empiricial cdfs and estimate zhat
# Later truncated observations, zero entries, will be indicated by x
eFx  <- apply(x,2,ecdf)
eFxx <- Map(function(f,x) do.call(f, list(x)), eFx, alply(x,2)  )
zhat <- matrix( unlist( lapply(eFxx,function(pr) qnorm( ( n/(n+1) )*pr ) ) ), n, p)


########################################################################
########   Initial values for gibbs sampling ###########################
########################################################################

		R_mc     <- cor(zhat) # Initial correlation matrix
	
		zhat_mc  <- zhat # Initial truncated observations with oberseved ones
	
		delta_mc <- qnorm( colSums( x==0 )/n ) # Initial values for Gaussian thresholds

		v0_mc   <- 0.001 # Initial value for spike variance

		tau_mc  <- h*v0_mc*matrix(1,p,p) # Initial spike variance matrix

		pijk_mc <- matrix( 2 / (p - 1),p,p) # Initial edge inclusion probabilities

        sig2_mc <- 5 # Initial tree scale

		U_mc <- rmvnorm(2, rep(0,p), sig2_mc*H_t ) # Initial latent positions

		
        burnin = 25000
        nmc = 100000

		bt <- Sys.time()
		set.seed(77843,kind = "Mersenne-Twister" ,sample.kind = "Rejection" )
		seed.state <- .Random.seed
		gibbsSample <- treeGibbs(x, delta_mc, zhat_mc, R_mc, v0_mc, tau_mc, pijk_mc, U_mc, sig2_mc, 
		                         hyperparameters, burnin, nmc, verbose=FALSE, sample.z=TRUE,thin=10)
		et <- Sys.time()
    print(et-bt)

    saveRDS(gibbsSample,file=paste(resultPath,"result.QMPchain_1.rds",sep=""))
print(.Random.seed)
