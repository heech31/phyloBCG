################################################################################ 
## Gibbs sampler for terminal nodes using Albert-Chip data augmentation       ##
## Author: Hee Cheol Chung                                                    ##
## Data: 08/08/2020                                                           ## 
################################################################################
library(truncnorm)
	# Uold <- U_mc
	# Vold <- V_mc
	# E <- E_mc



ACsample <- function(U,V,H,E, burnin, nmc){
	
	Uold <- U
	Vold <- V
	
	p <- ncol(Uold)
	K <- nrow(Uold)
	r <- dim(H)[1]
	
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
		
		ind_noj <- ind_noj_all[,jj]
		
		ejj  <- E[ind_noj,jj]
		
		uujj <- as.vector( crossprod(Uold[,-jj], Uold[,jj]) )
		

		yjj  <- rtruncnorm( p-1, a = 1-Inf^(1-ejj) , b = Inf^ejj - 1 , 
					mean = uujj, sd= 1 )

	    full_ind_noj <- c(ind_noj,(p+1):r)

		H22 <- H[full_ind_noj,full_ind_noj]

		H21 <- H[full_ind_noj, jj]
	
		uv <- cbind( Uold[,-jj], Vold[,])
	
		Phij   <- as.numeric( H[jj,jj] - t(H21) %*% solve(H22, H21) )
	
		thetaj <- as.vector( t(H21)%*%solve( H22 , t(uv) ) )
	
		Delj   <- solve( tcrossprod(Uold[,-jj],Uold[,-jj]) + diag(1/Phij,2) )

		gammaj <- Delj %*% ( Uold[,-jj]%*%yjj + diag(1/Phij,2)%*%thetaj )

		ujnew  <- rmvnorm(1,gammaj,Delj)

		Uold[,jj] <- ujnew

		}


    
    if (iter > burnin) {
      
      U_save[, , iter - burnin] = Uold
      
    }
}    
    
	results <- list(Unew = U_save)
 	return(results)
 	
 	
}

# Phi <- H %x% diag(K)

    # full_ind_noj <- c(ind_noj,(p+1):r)

	# H22 <- H[full_ind_noj,full_ind_noj]

	# H21 <- H[full_ind_noj, jj]
	
	# uv <- cbind( Uold[,-jj], Vold[,])
	
	# Phij   <- as.numeric( H[jj,jj] - t(H21) %*% solve(H22, H21) )
	
	# thetaj <- as.vector( t(H21)%*%solve( H22 , t(uv) ) )
	
	# Delj   <- solve( tcrossprod(Uold[,-jj],Uold[,-jj]) + diag(1/Phij,2) )

	# gammaj <- Delj %*% ( Uold[,-jj]%*%yjj + diag(1/Phij,2)%*%thetaj )

	# ujnew  <- rmvnorm(1,gammaj,Delj)

	# Uold[,jj] <- ujnew
	
	
	
      # Sig11 = Sig[ind_noi, ind_noi]
      # Sig12 = Sig[ind_noi, i]
      
      
      # invC11 = Sig11 - Sig12 %*% t(Sig12) / Sig[i, i]
      
      
      # Ci = (S[i, i] + lambda) * invC11 + diag(1. / tau_temp)
      
      
      # Ci = (Ci + t(Ci)) / 2
      
      # Ci_chol = chol(Ci)
      
      # mu_i = -solve(Ci_chol, solve(t(Ci_chol), S[ind_noi, i]))
      
      # beta = mu_i + solve(Ci_chol, rnorm(p - 1))
      
      
      
      # -solve(Ci, S[ind_noi, i])/ mu_i
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
 