rm(list=ls())
#setwd("/Users/heecheolchung/Dropbox/Research/TAMU/graphical/p50cdf/")
#funcPath   <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/functions/"
#resultPath <- "/Users/heecheolchung/Dropbox/Research/TAMU/graphical/p50cdf/"

setwd("/general/home/hcchung/graphical/p50cdf/m2/")
funcPath   <- "/general/home/hcchung/graphical/functions/"
resultPath <- "/general/home/hcchung/graphical/p50cdf/results/"

##Load required libraries
source(paste(funcPath,"libraries.R",sep=""))

ncores <- 10 # Number of cores used for parallel Gibbs sampling for 50 replicated datasets
registerDoParallel(ncores)

source(paste(funcPath,"genTree.R",sep=""))
source(paste(funcPath,"conditionalSample.R",sep=""))
source(paste(funcPath,"SSVS.R",sep=""))
source(paste(funcPath,"AlbertChib_exp.R",sep=""))
source(paste(funcPath,"distGibbs.R",sep=""))

nrep   <- 50  # Number of replicated datasets
burnin <- 500 # Number of burnin-iteration
nmc    <- 5000# Number of MCMC-iteration will be kept

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
#plot(graph_from_adjacency_matrix( Wtrue, mode="undirected"),layout = layout_with_kk)


########################################################################
########   Parameters for precision matrix sampling ####################
########################################################################
## Fix some hyperparameters
# Note that the parameterization used in the code is slightly different from those in Wang (2014).  )
h = 2500      # (v0 in code) = (v0 in paper)^2
lambda = 1    # Rate paramter of exponential prior for SSVS
IGsig2 = rep(1e-3,2) # Inv-gamma parameters for tree scale parameter (sigma2)
IGv0   = rep(1e-3,2) # Inv-gamma parameters for the spike and slab (v0)
hyperparameters <- list(h=h, lambda=lambda, IGsig2=IGsig2, IGv0=IGv0)


    



bt <- Sys.time()
########################################################################
########           Data Generation                  ####################
########################################################################

x50 <- vector("list",nrep)
for( ii in 1:nrep){
    x50[[ii]] <- synthData_from_ecdf(QMP, mar = 2, SigmaTrue , n=n, seed = NULL, verbose = FALSE)
  }
# Run MCMC for each 50 replicated data sets
Rep50 <- foreach( irep = 1:nrep)%dopar%{
  
    x <- x50[[irep]]


    eFx  <- apply(x,2,ecdf) # Obtain marginal cdf
    eFxx <- Map(function(f,x) do.call(f, list(x)), eFx, alply(x,2)  ) # Get F(x_j)
    zhat <- matrix( unlist( lapply(eFxx,function(pr) qnorm( ( n/(n+1) )*pr ) ) ), n, p) # Get \hat{z}_{o}


    ########################################################################
    ########   Initial values for gibbs sampling ###########################
    ########################################################################

    R_mc     <- cor(zhat)               # Initial correlation matrix
    zhat_mc  <- zhat                    # Initial truncated data. Observed data will be fixed
    delta_mc <- qnorm( colSums( x==0 )/n ) # Initial threshold
    v0_mc   <- 0.01                     # Initial spike variance
    tau_mc  <- h*v0_mc*matrix(1,p,p)    # Matrix representation of spike and slab variances
    pijk_mc <- matrix( 2 / (p - 1),p,p) # Initial edge inclusion probability
    djk  <- mytree$dist.terminal.nodes  # Tree distance of the terminal nodes
    gamma_mc <- 1                       # Scale of the distances


    gibbsSample <- distGibbs(x, delta_mc, zhat_mc, R_mc, v0_mc, tau_mc, pijk_mc, djk, gamma_mc,
                         hyperparameters, burnin, nmc, verbose=FALSE)


    Wtmp <- 1*( (Reduce("+" , alply(gibbsSample$E_gibbs[,,],3) )/dim(gibbsSample$E_gibbs)[3]) > 0.5 )

    edgeHat  <- Wtmp[ upper.tri(Wtmp) ]
    edgeTrue <- as.matrix(Wtrue)[ upper.tri(Wtmp) ] 
    edgeTotal <- sum( edgeTrue==1 )

    TP  <- sum( edgeHat[ edgeTrue==1 ] )
    FP  <- sum( edgeHat[ edgeTrue==0 ] )
    MCC <- cor(edgeHat,edgeTrue)


    OmegaHat <- ( Reduce("+" , alply(gibbsSample$C_gibbs[,,],3) )/nmc )
    What     <- ( Reduce("+" , alply(gibbsSample$E_gibbs[,,],3) )/nmc )
    pihat    <- ( Reduce("+" , alply(gibbsSample$pi_gibbs[,,],3) )/nmc )



    gibbsSummary <- list( 
                        What = What,
                        Wtrue = Wtrue,
                        edgeHat = edgeHat,
                        edgeTrue = edgeTrue,
                        edgeTotal = edgeTotal,
                        TP = TP,
                        FP = FP,
                        MCC = MCC,
                        OmegaHat = OmegaHat,
                        OmegaTrue = OmegaTrue,
                        pihat = pihat,
                        v0_gibbs = gibbsSample$v0_gibbs,
                        gamma = summary(gibbsSample$gamma_gibbs),
                        zsample.fail.count = gibbsSample$zsample.fail.count
                        )


    return( gibbsSummary )

    }
et <- Sys.time()

print(et-bt)

save.image(paste(resultPath,"Gibbs_cdf_D_",iseed,".RData",sep=""))




