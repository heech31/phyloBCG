################################################################################ 
## Gibbs sampler for the latent positions using Albert-Chip data augmentation ##
## Author: Hee Cheol Chung                                                    ##
## Date: 08/08/2020                                                           ##
## Update: 03/30/2020                                                         ##
################################################################################
	# Uold <- U_mc
	# Vold <- V_mc
	# E <- E_mc



ACsample_mrg <- function(U, H_t, E, burnin, nmc){
# U : K x p matrix with the latent positions in column
# H_t : tree covariance matrix without the scale parameter
# E   : Adjacency matrix
# burnin : Number of burnin iteration
# nmc    : Number of MCMC draws to keep
  
	Uold <- U #Latent positions  K x p 

	p <- ncol(Uold) # Number of variables
	K <- nrow(Uold) # The dimension of the latent space
	#r <- dim(H_t)[1] # Number of nodes. Should be the same as p
	
	U_save <- array(0, c(K, p, nmc))
	
	ind_noj_all <- matrix(0, p - 1, p)  
  
	for (jj in 1:p) {
		ind_noj_all[, jj] = setdiff(1:p,jj)
	}
	



for (iter in 1:(burnin + nmc)) {
  if (iter %% 2000 == 0) {
    print(paste("iter = ", iter, " nedge = ", (sum(Z) - p) / 2, sep = ""))
  }

	
	for( jj in 1:p ){
		
		ind_noj <- ind_noj_all[,jj] #Indices without j
		
		ejj  <- E[ind_noj,jj] # Edge indicators without jth diagonal component 
		# Uold is K x p 
		uujj <- as.vector( crossprod(Uold[,-jj, drop=FALSE], Uold[,jj, drop=FALSE]) ) # T_{-j}'*tj in the manucript
		
    # a and b are the lower and upper bound, respectively
		# If ejj = 0 then a = -Inf and b = 0 
		# If ejj = 1 then a = 0    and b = Inf 
		yjj  <- rtruncnorm( p-1, a = 1-Inf^(1-ejj) , b = Inf^ejj - 1 , 
					mean = uujj, sd= 1 )

	  full_ind_noj <- ind_noj

		H22 <- H_t[full_ind_noj,full_ind_noj] # Covariance matrix of the variables with the jth one.

		H21 <- H_t[full_ind_noj, jj] # Covariance between the jth variable and the rest.
	
		uv <- Uold[,-jj, drop=FALSE] # the latent position of the jth variable
	
		Phij   <- as.numeric( H_t[jj,jj] - t(H21) %*% solve(H22, H21) ) # Conditional variance of tj given T_{-j}
	
		thetaj <- as.vector( t(H21)%*%solve( H22 , t(uv) ) ) # Conditional mean of tj given T_{-j}
	  # Full conditional covariance matrix
		Delj   <- solve( tcrossprod(Uold[,-jj,drop=FALSE],Uold[,-jj,drop=FALSE]) + diag(1/Phij,K) ) 
		# K components of tj are independent fo Phi_j = diag(1/Phij,2)
    # Full condition mean vector
		gammaj <- Delj %*% ( Uold[,-jj]%*%yjj + diag(1/Phij,K)%*%thetaj )

		ujnew  <- mvtnorm::rmvnorm(1,as.vector(gammaj),Delj) # New jth latent position

		Uold[,jj] <- as.vector(ujnew) # Update the jth latent position
  
	}
    
  if (iter > burnin) {
    U_save[, , iter - burnin] = Uold
  }
  
}    
    
	results <- list(Unew = U_save)
 	return(results)
 	
 	
}

