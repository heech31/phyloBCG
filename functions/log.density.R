## This function calculates p(z_i| sth MCMC sample )
## Arithmetic mean is the posterior predictive ordinate, 
# PPO = mean_s{ p(z_i| sth MCMC sample ) }
## Harmonic mean is the conditional predictive ordiante 
# CPO = mean_s { p(z_i| sth MCMC sample )^{-1} }^{-1}
log.gaussian.density <- function(x, zhat_mc, delta_mc, R_mc){
  
  n <- nrow(x); p <- ncol(x);
  # Flag indicating observed variables of x
  o.flag.matrix <- x != 0
  log.density.i <- rep(0,n)
  n <- 10
  for( ii in 1:n){
    o.flag.i  <- o.flag.matrix[ii,]
    delta.o.i <- delta_mc[o.flag.i] # Lower threshold for observed variables of z[i,]
    # Partition correlation matrix
    R.oo <- R_mc[ o.flag.i, o.flag.i] # Variance-Covariance matrix of observed variables
    R.to <- R_mc[!o.flag.i, o.flag.i] # Covariance matrix between observed and truncated variables
    R.tt <- R_mc[!o.flag.i,!o.flag.i] # Variance-Covariance matrix of truncated variables
    
    # Conditional mean and variance of z.t | z.o, z.t < delta.t
    iR.tt.to <- solve(R.tt, R.to) 
    c.mu    <- as.vector( crossprod( zhat_mc[ii,!o.flag.i], iR.tt.to) )
    c.Sigma <- R.oo - crossprod(R.to, iR.tt.to)
    c.Sigma <- (c.Sigma + t(c.Sigma))/2
    
    # Sample z.t from truncated normal to complete z.new[-j]
    log.density.i[ii] <- dtmvnorm(zhat_mc[ii, o.flag.i], mean=c.mu, sigma= c.Sigma, 
                         lower= delta.o.i, upper = rep(Inf,length(c.mu)), log=TRUE   )
  }
  
  return(log.density.i)
  
}

