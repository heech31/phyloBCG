# sample from conditional distribution (truncated normal)
# This function samples truncated z conditioning on observed z.
zsample <- function(x,z,R,delta,below=TRUE){
	#library(tmvtnorm)
	# n     : Sample size
	# x     : p x 1 Observation vector with some zeros (not all)
	# z     : p x 1 Latent level normal vector ( f(z) = x )
	# R     : p x p Correlation matrix 
	# delta : p x 1 Truncation threshold for z (given).
    # below : Logical. If TRUE x<c for some c, and x>c, otherwise.

    ind0   <- (x==0) # Flag for 0 observations #sum(ind0)
	
	if( sum(ind0) == 0 ){
		warning("No truncation exists")
		results <- list(ind0 = ind0, z0new = NULL, cmu = NULL, cSigma = NULL)	
		return(results)
		break
	}
	
	
	delta0 <- delta[ind0] # truncation (censoring) thresholds corresponding to zero entries
	p0 <- sum(ind0) # The number of censored variables
	px <- sum(!ind0) # The number of observed (non-zero) variables
	z1 <- z[!ind0]   # The vector with observed z using observed x
	
	# Partition correlation matrix
	S00 <- matrix( R[ ind0, ind0], p0, p0) # Variance-Covariance matrix of truncated variables
	S01 <- matrix( R[ ind0,!ind0], p0, px) # Covariance matrix between truncated and truncated variables
	S11 <- matrix( R[!ind0,!ind0], px, px) # Variance-Covariance matrix of observed variables

	iS1101 <- solve(S11,t(S01) ) # invS11 * S10
	
	cmu    <- as.vector( crossprod( z1, iS1101  ) )  # Conditional mean S01*invS11*z1
	cSigma <- S00 - S01%*%iS1101           # Conditional covariance S00 - S01*invS11*S10
	cSigma <- ( cSigma + t(cSigma) )/2
	lowerb <- rep(-Inf,p0) # lower limits for truncated normal
	upperb <- delta0 # upper limits for truncated normal

	if(below){ # If truncated below
		upper <- delta0
		lower <- rep(-Inf,p0)
	}else{ # If truncated above
		upper <- rep(Inf,p0)
		lower <- delta0
	}

  zt.new  <- tmvtnorm::rtmvnorm(1, mean=cmu, sigma=cSigma, upper=upper, lower=lower,
	                        		 algorithm = "gibbs", burn.in.samples=50, thinning= 1 )

  # This is to check if conditional sampling fails or not
  # That happens when variance is too small and truncation point is far from the mean
  # For example, N(3,0.01), truncated above at -2, sampling from the truncated normal can be problematic
  # because the density below -2 is too thin, and rtmvnorm return NA.
  # The below part is added to monitor if that is happening and when it is happening.
  # After burn-in iteration, this should not happen.
  
  if( sum(is.na(zt.new))>0 ){ # If NA occurs
    zt.new  <- tmvtnorm::rtmvnorm(1, mean=cmu, sigma=(cSigma+diag(0.1,p0)), upper=upper, lower=lower,
                                algorithm = "gibbs", burn.in.samples=50, thinning= 1 )
    occur.fail <- 1
  }else{
    occur.fail <- 0
  }
	  

	results <- list(ind0 = ind0, zt.new = zt.new, occur.fail = occur.fail )
	return(results)
}




  
