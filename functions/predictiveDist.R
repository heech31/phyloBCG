# This function is for sampling from posterior predictive distribution xnew|xold
sample_post_predictive <- function(eFx, n, n.pred, delta_upper=NULL, R_mc){
  # eFx   : Scaled empirical cdf (list of length p)
  # n     : Sample size that is used for estimating eFx
  # delta_upper : Gaussian level threshold vector (not MCMC sample)
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

