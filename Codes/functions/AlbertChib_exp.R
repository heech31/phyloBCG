################################################################################ 
## Gibbs sampler for the tree scale parameter (gamma) of Dist model           ##
## Author: Hee Cheol Chung                                                    ##
## Data: 08/08/2020                                                           ## 
################################################################################

library(truncdist)#install.packages("truncdist")

# E <- E_mc
# gamma_mc <- 1
# gamma <- gamma_mc

ACsample_exp <- function(dij, E, gamma, burnin, nmc){
	
	utriInd <- upper.tri(dij)
	dij.vec <- dij[ 	utriInd ]
    eij.vec <- E[ 	utriInd ]
	
	p <- ncol(E)
	if( p*(p-1)/2 != length(dij.vec) ){
		prin("Error: dimensions do not match")
		break
	}
	
	pC2 <- p*(p-1)/2 
	
	gamma_save <- rep(0,nmc)
	
	lowerlim <- (0)^(1-eij.vec)
	
	upperlim <- Inf^eij.vec
	
  for (iter in 1:(burnin + nmc)) {

 
		# Sampling yij
		yij  <- rep(0,pC2)

		for( jj in 1:pC2){

			yij[jj]  <- rtrunc( n=1, spec="gamma", a = lowerlim[jj], b = upperlim[jj], shape=1, rate = gamma*dij.vec[jj] )

		}
		
		# Sampling gamma
		
		gamma <- rgamma(1, shape = pC2+1, rate = 1 + sum(dij.vec*yij) )
		
    
    if (iter > burnin) {
      
      gamma_save[iter - burnin] = gamma
      
    }
}    
    
	results <- list(gamma_save = gamma_save)
 	return(results)
 	
 	
}



