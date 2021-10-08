
quantile.stepcdf <- function(prob, empf, interval, tol = 1e-3, maxiter = 100){
  ans <- c()
  # uniroot.all from rootSolve package was the fastest one.
  sol <- as.numeric(
          rootSolve::uniroot.all(function(x){empf(x)-prob}, interval = interval, tol = tol, maxiter = maxiter))
  
  if( length(sol) != 0 ){
    # in case quantile function is a step function, added the following step.
    if (prob <= empf(floor(sol))) {
      ans <- floor(sol)
    } else if (prob > empf(floor(sol))) {
      ans <- ceiling(sol)
    }
  }else{
    ans <- 0
  }
  return(ans)
}







conditional_predictive_mean <- function(eFx, x.new, delta_mc, R_mc, maxbound){
  # eFx   : Scaled empirical cdf (list of length p)
  # x.new : Test observation that is not used for model fitting (vector of length p)
  # delta_mc : Current MCMC sample of threshold vector
  # R_mc     : Current MCMC sample of correlation matrix
  
  # Flag indicating observed variables of x.new
  o.flag <- x.new!=0
  # Probability of Xj less than equal to xj
  F.new  <- unlist( Map(function(f,x) do.call(f, list(x)), eFx, alply(x.new,1) ) )
  # If prob is 0 then replace it to 0.01 
  F.new[ F.new==0 ] <- 0.01
  # Estimated Gaussian vector z.new
  z.new  <- qnorm(F.new)
  z.new[!o.flag] <- 0 #
  
  # z.new[j] given z.new[-j] and x will be stored for j=1,...,p
  # x.new[j] given z.new[-j] and x will be stored for j=1,...,p
  z.new.loo <- rep(0,p)
  z.new.loo.trc <- rep(0,p)
  x.new.loo <- rep(0,p)
  
  for( jj in 1:p){
    
    R.no.j       <- R_mc[-jj,-jj] # Correlation matrix without the jth variable
    z.new.no.j   <- z.new[-jj]    # z.new[-j]
    o.flag.no.j  <- o.flag[-jj]   # Flag for non-zeros variables of z.new[-j]
    delta.t.no.j <- delta_mc[-jj][!o.flag.no.j] # Threshold for truncated variables of z.new[-j]
    # Partition correlation matrix
    R.oo <- R.no.j[ o.flag.no.j, o.flag.no.j] # Variance-Covariance matrix of observed variables
    R.ot <- R.no.j[ o.flag.no.j,!o.flag.no.j] # Covariance matrix between observed and truncated variables
    R.tt <- R.no.j[!o.flag.no.j,!o.flag.no.j] # Variance-Covariance matrix of truncated variables
    
    # Conditional mean and variance of z.t | z.o, z.t < delta.t
    iR.oo.ot <- solve(R.oo, R.ot) 
    c.mu    <- as.vector( crossprod( z.new.no.j[o.flag.no.j], iR.oo.ot) )
    c.Sigma <- R.tt - crossprod(R.ot, iR.oo.ot)
    c.Sigma <- (c.Sigma + t(c.Sigma))/2
    
    # Sample z.t from truncated normal to complete z.new[-j]
    from.truncated.mvnorm <- rtmvnorm(1, mean=c.mu, sigma= c.Sigma, 
                                      lower= rep(-Inf,length(c.mu)), upper = delta.t.no.j,
                                      algorithm="gibbs",burn.in.sample=10)
    
    z.new.no.j[!o.flag.no.j] <- as.vector(from.truncated.mvnorm)
    
    # Conditional mean and standard deviation of z.new[j] | z.new[-j]
    mu.jj.given._jj <- sum( solve( R.no.j,R_mc[jj,-jj] )*z.new.no.j )
    sd.jj.given._jj <- sqrt( R_mc[jj,jj] - sum( solve( R.no.j, R_mc[jj,-jj] ) * R_mc[jj,-jj] ) )
    
    z.new.loo[jj] <- rnorm(1,mu.jj.given._jj, sd.jj.given._jj)
    z.new.loo.trc[jj] <- ifelse( z.new.loo[jj]<delta_mc[jj], 0, z.new.loo[jj] )
    
    if(z.new.loo[jj] != 0){
      prob <- pnorm( z.new.loo[jj] )
      interval <- c(0,maxbound[jj] )
      x.new.loo[jj]  <- quantile.stepcdf(prob, eFx[[jj]], interval )
      }
    
  }#End jj
  
  res <- list(x.new.loo = x.new.loo, z.new.loo = z.new.loo, z.new.loo.trc = z.new.loo.trc)
  
  return(res)
  
}







