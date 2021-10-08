
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







sample_post_predictive <- function(eFx, n, n.pred, delta_upper=NULL, R_mc){
  # eFx   : Scaled empirical cdf (list of length p)
  # n     : Sample size that is used for estimating eFx
  # delta_mc : Current MCMC sample of threshold vector
  # R_mc     : Current MCMC sample of correlation matrix

  z.pred  <- mvtnorm::rmvnorm(n.pred, mean=rep(0,ncol(R_mc)), sigma = R_mc )
  x.pred  <- z.pred*0

  for( jj in 1:ncol(R_mc) ){
    # Check if the value is less than the corrosponding threshold
    nonzero.ind <- z.pred[,jj] > delta_upper[jj]
    # Get the probability (normal cdf)
    prob <- pnorm( z.pred[nonzero.ind,jj] )
    # Interval where the solution will be searched
    interval <- c(0,maxbound[jj] )
    # Vectorize pseudo inverse of emprical cdfs
    qs.vec <- Vectorize(function(probs) quantile.stepcdf(probs, eFx[[jj]], interval ) )
    # Apply pseudo inverse of emprical cdfs
    if(n.pred>1){
      x.pred[nonzero.ind,jj]  <- qs.vec(prob)
    }else{
      if(nonzero.ind){
        x.pred[1,jj]  <- qs.vec(prob)
      }
    }
  }

  samples <- list(z.pred = z.pred, x.pred = x.pred)

  return(samples)
}


# 
# 
# which.one.is.close <- function(x,y){
#   min.ind <- which.min( abs(x-y) )
#   return(min.ind)
# }
# 
# 
# 
# 
# sample_post_predictive <- function(eFx, n.pred, delta_upper=NULL, R_mc){
#   # eFx   : Scaled empirical cdf (list of length p)
#   # delta_mc : Current MCMC sample of threshold vector
#   # R_mc     : Current MCMC sample of correlation matrix
#   
#   z.pred  <- mvtnorm::rmvnorm(n.pred, mean=rep(0,ncol(R_mc)), sigma = R_mc )
#   x.pred  <- z.pred*0
#   
#   for( jj in 1:ncol(R_mc) ){
#     nonzero.ind <- z.pred[,jj] > delta_upper[jj]
#     prob <- ((n.pred+1)/n.pred)*pnorm( z.pred[nonzero.ind,jj] )
#     xseq <- seq(0,maxbound[jj],l=1000)
#     smoothed.cdf <- smooth.spline(xseq, eFx[[jj]](xseq), df=20)
#     
#     x.ind <- Vectorize( function(x) which.one.is.close( smoothed.cdf$y,x) )(prob )
#     x.pred[nonzero.ind,jj]  <- xseq[x.ind]
#   }
#   
#   samples <- list(z.pred = z.pred, x.pred = x.pred)
#   
#   return(samples)
# }
# 

