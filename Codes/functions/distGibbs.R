##########################################################
## Gibbs sampler for the tree distance model            ##
## Author: Hee Cheol Chung                              ##
## Data: 03/08/2020                                     ##
##########################################################

distGibbs <- function(x, delta_mc, zhat_mc, R_mc, v0_mc, tau_mc, pijk_mc, djk, gamma_mc,
                      hyperparameters, burnin, nmc, verbose=FALSE){
## x: zero-inflated observed data (n x p)
## delta_mc : Initial value of Gaussian thresholds (p x 1)
## zhat_mc  : Initial value of latent Gaussian data (n x p)
## R_mc     : Initial value of latent Gaussian correlation matrix (p x p)
## tau_mc   : Matrix with spike and slab variance componets (p x p)
## pijk_mc  : Initial value of edge inclusion probability matrix (p x p)
## djk      : Tree distances matrix (p x p) between node
## hyperparameters : List of hyperparameter values (h,lambda,IGsig2, IGv0)

h      <- hyperparameters$h
lambda <- hyperparameters$lambda
IGsig2 <- hyperparameters$IGsig2
IGv0   <- hyperparameters$IGv0

np <- dim(zhat_mc) # Sample size and dimension
n <- np[1]; p <- np[2]

# Gibbs sample of  
R_gibbs     <- array(0,dim=c(p,p,nmc)) # Correlation matrix
C_gibbs     <- array(0,dim=c(p,p,nmc)) # Concentration matrix
E_gibbs     <- array(0,dim=c(p,p,nmc)) # Adjacency matrix
delta_gibbs <- array(0,dim=c(p,1,nmc)) # Gaussian thresholds

gamma_gibbs <- rep(0,nmc)
v0_gibbs    <- rep(0,nmc)
tau_update  <- array(0,dim=c(p,p,nmc)) 
pi_update   <- array(0,dim=c(p,p,nmc))

utriFlag <- upper.tri(matrix(NA,p,p)) ## Upper triangular flag
zsample.fail.count <- integer(length = nmc)

for( gg in 1:(burnin+nmc)){

########################
#### Sample zhat_0  ####
########################
#cbind(x[ii,],zhat_mc[ii,],delta_mc)
fail.count <- 0
for( ii in 1:n){# Sample truncated Gaussian variables
tmp  <- zsample(x[ii,],zhat_mc[ii,],R_mc, delta_mc)
zhat_mc[ii,tmp$ind0] <- as.vector( tmp$zt.new )
if( sum(is.na(tmp$zt.new))>0 ){
    print("NA occured during zsampling");
    break} 
    fail.count <- fail.count + tmp$occur.fail
}
if(gg>burnin){
zsample.fail.count[gg] <- as.integer(fail.count)
}
########################
#### Sample delta   ####
########################

zmax_mc <- rep(0,p) # variables-wise maximum among the truncated observations (sampled)
zmin_mc <- rep(0,p) # variables-wise minimum among the non-zero observations (observed)
for( jj in 1:p){
  zmax_mc[jj] <- max( zhat_mc[x[,jj]==0,jj] )# 
  zmin_mc[jj] <- min( zhat_mc[x[,jj]!=0,jj] )
}

if( sum(zmax_mc>zmin_mc)!=0  ){ warning("zmin has to be larger than zmax")}#cbind(zmax_mc,zmin_mc)
delta_mc <- runif(p,zmax_mc,zmin_mc) # Sample thresholds from uniform
delta_mc[is.nan(delta_mc)] <- -Inf   # For variable with no truncated obs

# zmax_mc = +inf gives NaN
# so make the threshold -Inf (no trunction)
########################################
#### Sample Sigma (Omega) and Graph ####
########################################

v1_mc = h * v0_mc # Slab variance

V0 = v0_mc * matrix(1, p, p); 
V1 = v1_mc * matrix(1, p, p);

S    <- crossprod(zhat_mc) # Gram matrix

SSSL <- SSVS(S, n, Sig = R_mc, V0, V1, tau_mc, lambda, pijk_mc) # One Gibbs iteration



R_mc <- cov2cor( SSSL$Sig_new )
C_mc <- SSSL$C_new
E_mc <- SSSL$Z_new
tau_mc <- SSSL$T_new


########################################
####   Sample spike varinace v0^2   ####
########################################


pp <- p*(p-1)/2
wijsq <- C_mc[utriFlag]^2

v0shape <- pp/2 + IGv0[1] # pp/2 = p(p-1)/4
v0rate  <- sum( wijsq/(h^(E_mc[utriFlag]) ) )/2 + IGv0[2]

v0_mc   <- 1/rgamma(1, shape=v0shape, rate=v0rate)



#################################################
#### Sample gamma ( pi_ij = exp(-gamma*djk) ) ###
#################################################

gamma_new <- ACsample_exp(djk, E_mc, gamma_mc, burnin = 0, nmc = 1)$gamma_save

gamma_mc <- gamma_new

pijk_mc  <- exp( - gamma_mc * djk )


















if( gg > burnin){
delta_gibbs[,,gg-burnin] <- delta_mc #
R_gibbs[,,gg-burnin] <- R_mc
C_gibbs[,,gg-burnin] <- C_mc
E_gibbs[,,gg-burnin] <- E_mc
v0_gibbs[gg-burnin]  <- v0_mc
tau_update[,,gg-burnin] <- tau_mc
pi_update[,,gg-burnin] <- pijk_mc

gamma_gibbs[gg-burnin] <- gamma_mc


}

if( verbose ){
if( gg%%2000 == 0 ){ print(gg);  print( Sys.time() )  }
}

}


gibbsSample <- list(
R_gibbs = R_gibbs,
C_gibbs = C_gibbs,
E_gibbs = E_gibbs,
pi_gibbs = pi_update,
v0_gibbs = v0_gibbs,
gamma_gibbs = gamma_gibbs,
zsample.fail.count = zsample.fail.count
)



}




