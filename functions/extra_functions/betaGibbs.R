##############################################
## Gibbs sampler for the beta model         ##
## Author: Hee Cheol Chung                  ##
## Data: 03/08/2020                         ##
##############################################

betaGibbs <- function(x, x.new, delta_mc, zhat_mc, R_mc, v0_mc, h, tau_mc, pi_mc,
                      hyperparameters, burnin, nmc, verbose=FALSE){
## x: zero-inflated observed data (n x p)
## delta_mc : Initial value of Gaussian thresholds (p x 1)
## zhat_mc  : Initial value of latent Gaussian data (n x p)
## R_mc     : Initial value of latent Gaussian correlation matrix (p x p)
## tau_mc   : Matrix with spike and slab variance componets (p x p)
## pi_mc    : Initial value of the edge inclusion probability (0,1)
## hyperparameters : List of hyperparameter values (h,lambda,IGsig2, IGv0)



h      <- hyperparameters$h
lambda <- hyperparameters$lambda
IGsig2 <- hyperparameters$IGsig2
IGv0   <- hyperparameters$IGv0

np <- dim(zhat_mc) # Sample size and dimension
n <- np[1]; p <- np[2]







## This part is added 08/04/2021 by Hee Cheol Chung
###########for predictive means ################
if( !is.null(x.new) ){
  ecdf.scale <- n/(n+1)
  eFx  <- apply(x,2,ecdf)
  eFx  <- lapply(eFx, function(x){  function(y)  ecdf.scale *x(y) })
  maxbound <- apply(x, 2, max) # to restrict the search range of the solution (for uniroot.all)
  x.loo_gibbs <- matrix(0,p,nmc)
  z.loo_gibbs <- matrix(0,p,nmc)
}else{
  x.loo_gibbs <- matrix(NA,1,nmc)
  z.loo_gibbs <- matrix(NA,1,nmc)
}
################################################


# Gibbs sample of
R_gibbs     <- array(0,dim=c(p,p,nmc)) # Correlation matrix
C_gibbs     <- array(0,dim=c(p,p,nmc)) # Concentration matrix
E_gibbs     <- array(0,dim=c(p,p,nmc)) # Adjacency matrix
delta_gibbs <- array(0,dim=c(p,1,nmc)) # Gaussian thresholds
pi_gibbs    <- rep(0,nmc)  # Edge includsion probability

v0_gibbs    <- rep(0,nmc) # Spike variance
tau_update  <- array(0,dim=c(p,p,nmc)) # Spike and slab variance entries
pi_update   <- array(0,dim=c(p,p,nmc)) # Edge inclusion probabilty matrix
z_gibbs     <- array(0,dim=c(n,p,nmc)) # Latent positions
#PDM_gibbs   <- matrix(0, p, nmc) # Pivotal discrepancy measure
#log.den_gibbs <- matrix(0,n,nmc) # log.gaussian.densities
utriFlag <- upper.tri(matrix(NA,p,p)) ## Upper triangular flag
zsample.fail.count <- integer(length = nmc)

for( gg in 1:(burnin+nmc)){

########################
#### Sample zhat_0  ####
########################
  
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
  zmax_mc[jj] <- max( zhat_mc[x[,jj]==0,jj] )
  zmin_mc[jj] <- min( zhat_mc[x[,jj]!=0,jj] )
}

if( sum(zmax_mc>zmin_mc)!=0  ){ warning("zmin has to be larger than zmax")}
delta_mc <- runif(p,zmax_mc,zmin_mc) # Sample thresholds from uniform
delta_mc[is.nan(delta_mc)] <- -Inf   # For variable with no truncated obs
# zmax_mc = +inf give NaN
# so make the threshold -Inf (no trunction)
########################################
#### Sample Sigma (Omega) and Graph ####
########################################

V0 = v0_mc * matrix(1, p, p)
V1 = h * v0_mc * matrix(1, p, p)

S    <- crossprod(zhat_mc) # Gram matrix

SSSL <- SSVS(S, n, Sig = R_mc, V0, V1, tau_mc, lambda, pi_mc) # One Gibbs iteration



R_mc <- cov2cor( SSSL$Sig_new )
R_mc <- ( R_mc + t(R_mc) )/2
C_mc <- SSSL$C_new
E_mc <- SSSL$Z_new



########################################
####   Sample spike varinace v0^2   ####
########################################


pp <- p*(p-1)/2
wijsq <- C_mc[utriFlag]^2

v0shape <- pp/2 + IGv0[1] # pp/2 = p(p-1)/4
v0rate  <- sum( wijsq/(h^(E_mc[utriFlag]) ) )/2 + IGv0[2]

v0_mc   <- 1/rgamma(1, shape=v0shape, rate=v0rate)

v1_mc = h * v0_mc # Slab variance
tau_mc[E_mc==0] <- v0_mc
tau_mc[E_mc==1] <- v1_mc


########################################
####         Sample pijk_mc         ####
########################################

ne1 <- sum(E_mc[utriFlag])
pi_mc <- rbeta(1,ne1,pp-ne1)



################################################################################################
############################################
#### Get x.new[j]|x.new[-j],x           ####
############################################
if( !is.null(x.new) ){
  predict.z.x <- conditional_predictive_mean(eFx, as.vector(x.new), delta_mc, R_mc, maxbound)
  x.new.loo_mc <- predict.z.x$x.new.loo
  z.new.loo_mc <- predict.z.x$z.new.loo
}else{
  x.new.loo_mc <- NA
  z.new.loo_mc <- NA
  }
################################################################################################


################################################################################################
############################################################
#### Get p(z_i| sth MCMC sample)  (for PPO and CPO)     ####
############################################################
#log.den_mc <- log.gaussian.density(x,zhat_mc,delta_mc,R_mc)
################################################################################################
                
        
if( gg > burnin){
delta_gibbs[,,gg-burnin] <- delta_mc
R_gibbs[,,gg-burnin] <- R_mc # Save correlation matrix
C_gibbs[,,gg-burnin] <- C_mc # Save concentration matrix
E_gibbs[,,gg-burnin] <- E_mc # Save ajacency matrix
v0_gibbs[gg-burnin]  <- v0_mc # Save spike variance
pi_gibbs[gg-burnin]  <- pi_mc
tau_update[,,gg-burnin] <- tau_mc # Updated spike and slab entries
z_gibbs[,,gg-burnin]  <- zhat_mc
x.loo_gibbs[,gg-burnin] <- x.new.loo_mc
z.loo_gibbs[,gg-burnin] <- z.new.loo_mc
#PDM_gibbs[,gg-burnin] <-
#  PDM( x = x, z = zhat_mc, R = R_mc, delta = delta_mc)
#log.den_gibbs[,gg-burnin] <- log.den_mc
}

if( verbose ){
if( gg%%2000 == 0 ){ print(gg);  print( Sys.time() )  }
}

}


gibbsSample <- list(
delta_gibbs = delta_gibbs,
R_gibbs = R_gibbs,
C_gibbs = C_gibbs,
E_gibbs = E_gibbs,
v0_gibbs = v0_gibbs,
pi_gibbs = pi_gibbs,
x.loo_gibbs = x.loo_gibbs,
z.loo_gibbs = z.loo_gibbs,
z_gibbs = z_gibbs,
#PDM_gibbs = PDM_gibbs,
#log.den_gibbs = log.den_gibbs,
zsample.fail.count = zsample.fail.count
)




}

