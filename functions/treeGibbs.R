##########################################################
## Gibbs sampler for the phylogenetic tree model        ##
## Author: Hee Cheol Chung                              ##
## Date: 03/08/2021                                     ##
##########################################################

treeGibbs <- function(x, delta_mc, zhat_mc, R_mc, v0_mc, tau_mc, pijk_mc, U_mc, sig2_mc,
                      hyperparameters, burnin, nmc, x.new=NULL, verbose=FALSE, thin=NULL){
## x: zero-inflated observed data (n x p)
## delta_mc : Initial value of Gaussian thresholds (p x 1)
## zhat_mc  : Initial value of latent Gaussian data (n x p)
## R_mc     : Initial value of latent Gaussian correlation matrix (p x p)
## tau_mc   : Matrix with spike and slab variance componets (p x p)
## pijk_mc  : Initial value of edge inclusion probability matrix (p x p)
## U_mc     : Initial value of latent positions (L x p)
## hyperparameters : List of hyperparameter values (h,lambda,IGsig2, IGv0)
## burnin   : Number of burn-in iteration
## nmc      : Number of MCMC ieration

h      <- hyperparameters$h
lambda <- hyperparameters$lambda
IGsig2 <- hyperparameters$IGsig2
IGv0   <- hyperparameters$IGv0
H_t    <- hyperparameters$H_t
H_tinv <- hyperparameters$H_tinv

np <- dim(zhat_mc) # Sample size and dimension
n <- np[1]; p <- np[2] 

if(is.null(thin)){
  nnmc <- nmc
}else{
  nnmc <- round(nmc/thin)
}

## This part is added 08/04/2021 by Hee Cheol Chung
###########for predictive means ################
if( !is.null(x.new) ){
  ecdf.scale <- n/(n+1)
  eFx  <- apply(x,2,ecdf)
  eFx  <- lapply(eFx, function(x){  function(y)  ecdf.scale *x(y) })
  maxbound <- apply(x, 2, max) # to restrict the search range of the solution (for uniroot.all)
  x.loo_gibbs <- matrix(0,p,nnmc)
  z.loo_gibbs <- matrix(0,p,nnmc)
}else{
  x.loo_gibbs <- matrix(NA,1,nnmc)
  z.loo_gibbs <- matrix(NA,1,nnmc)
}
################################################


# Gibbs sample of
R_gibbs     <- array(0,dim=c(p,p,nnmc)) # Correlation matrix
C_gibbs     <- array(0,dim=c(p,p,nnmc)) # Concentration matrix
E_gibbs     <- array(0,dim=c(p,p,nnmc)) # Adjacency matrix
delta_gibbs <- array(0,dim=c(p,1,nnmc)) # Gaussian thresholds
U_gibbs     <- array(0,dim=c(K,p,nnmc)) # Latent positions
sig2_gibbs  <- rep(0,nnmc) # Tree scale parameter
v0_gibbs    <- rep(0,nnmc) # Spike variance
tau_update  <- array(0,dim=c(p,p,nnmc)) # Spike or slab variances for omega_jk
pi_update   <- array(0,dim=c(p,p,nnmc)) # Edge inclusion probability for e_jk
z_gibbs     <- array(0,dim=c(n,p,nnmc)) # Latent positions
utriFlag <- upper.tri(matrix(NA,p,p))   # Upper triangular flag
zsample.fail.count <- integer(length = nnmc)
count <- 1 # For thinning index





for( gg in 1:(burnin+nmc)){## gg<-1

########################
#### Sample zhat_0  ####
########################
#gg=1  
fail.count <- 0
for( ii in 1:n){# Sample truncated Gaussian variables
tmp  <- zsample(x[ii,],zhat_mc[ii,],R_mc, delta_mc)
zhat_mc[ii,tmp$ind0] <- as.vector( tmp$zt.new )
if( sum(is.na(tmp$zt.new))>0 ){
    print("NA occured during zsampling");
    break} 
    fail.count <- fail.count + tmp$occur.fail
}



########################
#### Sample delta   ####
########################

zmax_mc <- rep(0,p) # variables-wise maximum among the truncated observations (sampled)
zmin_mc <- rep(0,p) # variables-wise minimum among the non-zero observations (observed)
for( jj in 1:p){
  zmax_mc[jj] <- max( zhat_mc[x[,jj]==0,jj] )
  zmin_mc[jj] <- min( zhat_mc[x[,jj]!=0,jj] )
}

if( sum(zmax_mc>zmin_mc)!=0  ){ warning("zmin has to be larger than zmax");}
delta_mc <- runif(p,zmax_mc,zmin_mc) # Sample thresholds from uniform
delta_mc[is.nan(delta_mc)] <- -Inf   # For variable with no truncated obs
# zmax_mc = +inf give NaN
# so make the threshold -Inf (no trunction)
########################################
#### Sample Sigma (Omega) and Graph ####
########################################

V0 = v0_mc * matrix(1, p, p) # v0*1_{p,p}
V1 = h * v0_mc * matrix(1, p, p) # v1*1_{p,p}

S    <- crossprod(zhat_mc) # Gram matrix

SSSL <- SSVS(S, n, Sig = R_mc, V0, V1, tau_mc, lambda, pijk_mc) # One Gibbs iteration



R_mc <- cov2cor( SSSL$Sig_new )
R_mc <- ( R_mc + t(R_mc) )/2
C_mc <- SSSL$C_new
E_mc <- SSSL$Z_new



########################################
####   Sample spike varinace v0^2   ####
########################################


pp <- p*(p-1)/2
wijsq <- C_mc[utriFlag]^2

v0shape <- pp/2 + IGv0[1]
v0rate  <- sum( wijsq/(h^(E_mc[utriFlag]) ) )/2 + IGv0[2]

v0_mc   <- 1/rgamma(1, shape=v0shape, rate=v0rate)

v1_mc = h * v0_mc # Slab variance
tau_mc[E_mc==0] <- v0_mc
tau_mc[E_mc==1] <- v1_mc


########################################
#### Sample latent positions        ####
########################################

# We assume that H_t is given ( fixed )
# Sample U_mc (K x p matrix of latent positions)
# Albert-Chip data augmentation. _mrg indicates that 
# we are only considering marginal covariance matrix of the terminal nodes.
# There was a version sampling internal and terminal nodes together.
U_mc <- ACsample_mrg(U_mc, sig2_mc*H_t, E_mc, burnin = 0, nmc = 1)$Unew[,,1]
U_mc <- matrix(U_mc,K,p)

Gshape <- K*p/2 + IGsig2[1] # Shape of IG gamma for tree scale
Grate  <- sum( diag( U_mc %*% tcrossprod(H_tinv, U_mc) ) )/2 + IGsig2[2] # Rate of IG gamma for tree scale
sig2_mc <- 1/rgamma(1, shape = Gshape,  rate = Grate)


############################################
#### Get piij (edge inclusion prob)     ####
############################################

simil           <- crossprod(U_mc)

# Update edge probabilities
pijk_mc <- pnorm(simil)

####################################################################
####################################################################
#### Sampling from posterior condtional predictive distribution ####
#### Get x.new[j]|x.new[-j],x                                   ####
####################################################################
if( !is.null(x.new) ){
  predict.z.x <- conditional_predictive_mean(eFx, as.vector(x.new), delta_mc, R_mc, maxbound)
  x.new.loo_mc <- predict.z.x$x.new.loo
  z.new.loo_mc <- predict.z.x$z.new.loo
}else{
  x.new.loo_mc <- NA
  z.new.loo_mc <- NA
  }
####################################################################




if( gg > burnin){
  if( is.null(thin) ){
      delta_gibbs[,,gg-burnin]<- delta_mc
      R_gibbs[,,gg-burnin]    <- R_mc
      C_gibbs[,,gg-burnin]    <- C_mc
      E_gibbs[,,gg-burnin]    <- E_mc
      v0_gibbs[gg-burnin]     <- v0_mc
      tau_update[,,gg-burnin] <- tau_mc
      pi_update[,,gg-burnin]  <- pijk_mc
      U_gibbs[,,gg-burnin]  <- U_mc
      sig2_gibbs[gg-burnin] <- sig2_mc
      z_gibbs[,,gg-burnin]  <- zhat_mc
      x.loo_gibbs[,gg-burnin] <- x.new.loo_mc
      z.loo_gibbs[,gg-burnin] <- z.new.loo_mc
      zsample.fail.count[gg-burnin] <- as.integer(fail.count)
  }else{
  if( (gg-burnin)%%thin == 0){
    delta_gibbs[,,count]<- delta_mc
    R_gibbs[,,count]    <- R_mc
    C_gibbs[,,count]    <- C_mc
    E_gibbs[,,count]    <- E_mc
    v0_gibbs[count]     <- v0_mc
    tau_update[,,count] <- tau_mc
    pi_update[,,count]  <- pijk_mc
    U_gibbs[,,count]  <- U_mc
    sig2_gibbs[count] <- sig2_mc
    z_gibbs[,,count]  <- zhat_mc
    x.loo_gibbs[,count] <- x.new.loo_mc
    z.loo_gibbs[,count] <- z.new.loo_mc
    zsample.fail.count[count] <- as.integer(fail.count)
    count <- count + 1
    }
  }
}

  if( verbose & gg%%2000 == 0 ){ print(gg);  print( Sys.time() )  }

}


gibbsSample <- list(
delta_gibbs = delta_gibbs,
R_gibbs = R_gibbs,
C_gibbs = C_gibbs,
E_gibbs = E_gibbs,
pi_gibbs = pi_update,
v0_gibbs = v0_gibbs,
U_gibbs = U_gibbs,
sig2_gibbs = sig2_gibbs,
x.loo_gibbs = x.loo_gibbs,
z.loo_gibbs = z.loo_gibbs,
z_gibbs = z_gibbs,
zsample.fail.count = zsample.fail.count
)

return(gibbsSample)
}





