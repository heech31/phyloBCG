PhyloBCG
================

This folder contains R codes that produce the simulation and real data
analysis results presented in the paper “Phylogenetically informed
Bayesian truncated copula graphical models for microbial association
networks” by Hee Cheol Chung, Irina Gaynanova and Yang Ni. [arXiv
link](https://arxiv.org/pdf/2105.05082.pdf)

# Simple example

### Load data, packages, and functions

``` r
sourcePath <- "./functions/"
dataPath <- "./QMPdata/"

# Load packages and source functions
files.sources <-  list.files(path=sourcePath)
files.sources <- setdiff( files.sources, "TreeDataGeneration.R" )
extension <- substring(files.sources,nchar(files.sources)-1,nchar(files.sources))
lets.source <- paste(sourcePath, files.sources[extension==".R"], sep="")
mapply(source, lets.source) 


#Load QMP data (need to be loaded after packages and functions)
load(paste(dataPath,"QMPtree.RData",sep="") )
```

### Summary of QMP data

``` r
# Summary
t( apply(QMP,2,summary) )

# Sample size and dimension
(n <- dim(QMP)[1])  # Sample size
(p <- dim(QMP)[2])  # Dimension
x <- QMP

# Zero-proportions
colMeans(QMP==0)

# Phylogenetic tree
plot(qmptree, type="phylogram")
nodelabels(qmptree$node.label, cex=0.8, adj=1)
```

![](./README_QMP_files/figure-gfm/example-1.png)<!-- -->

``` r
# Tree correlation matrix H
corrplot(H)
```

![](README_QMP_files/figure-gfm/example-2.png)<!-- -->

### Set hyperparameters

``` r
########################################################################
#####             Set hyperarameters                ####################
########################################################################
## Fix some hyperparameters
# Note that the parameterization used in the code is slightly different from those in Wang (2014).  )
h <-  2500      # (v0 in code) = (v0 in paper)^2
lambda <-  1    # Rate paramter of exponential prior for SSVS
IGsig2 <-  rep(1e-3,2) # Inv-gamma parameters for tree scale parameter (sigma2)
IGv0   <- rep(1e-3,2) # Inv-gamma parameters for the spike and slab (v0)
hyperparameters <- list(h=h, lambda=lambda, IGsig2=IGsig2, IGv0=IGv0)
K <- 2 # Latent space dimension
```

### Set initial values

``` r
# Get scaled empirical cdfs and zhat
ecdf.scale <- n/(n+1)
eFx  <- apply(x,2,ecdf)
eFx  <- lapply(eFx, function(x){  function(y)  ecdf.scale *x(y) })
eFxx <- Map(function(f,x) do.call(f, list(x)), eFx, alply(x,2)  )
zhat <- matrix( unlist( lapply(eFxx,function(pr) qnorm( pr ) ) ), n, p)
H_t <- cov2cor(H)
H_tinv <- solve(H_t)

########################################################################
########   Initial values for gibbs sampling ###########################
########################################################################

R_mc     <- cor(zhat)               # Initial correlation matrix
zhat_mc  <- zhat                    # Initial truncated data. Observed data will be fixed
delta_mc <- qnorm( colSums( x==0 )/n ) # Initial threshold
v0_mc   <- 0.01                     # Initial spike variance
tau_mc  <- h*v0_mc*matrix(1,p,p)    # Matrix representation of spike and slab variances
pijk_mc <- matrix( 2 / (p - 1),p,p) # Initial edge inclusion probability
U_mc <- mvtnorm::rmvnorm(K, rep(0,p), H_t )  # Initial latent positions
sig2_mc <- 1                        # Initial tree scale parameter
```

### Run treeGibbs

``` r
burnin <- 50  # Number of burnin-iteration
nmc    <- 50  # Number of MCMC-sample that will be kept

gibbsSample <- treeGibbs(x, x.new=NULL, delta_mc, zhat_mc, R_mc, v0_mc, tau_mc, pijk_mc, U_mc, sig2_mc,
                         hyperparameters, burnin, nmc, verbose=FALSE, thin=NULL)
```