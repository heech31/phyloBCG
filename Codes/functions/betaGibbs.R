
########################################################################
########################################################################
########              Gibbs sampling               #####################
########################################################################
########################################################################

betaGibbs <- function(x, zhat, delta_mc, zhat_mc, R_mc, v0_mc, h,
                      tau_mc, pijk, burnin, nmc, verbose){

	p <- dim(zhat)[2]; pp <- p*(p-1)/2

	R_gibbs     <- array(0,dim=c(p,p,nmc))
	C_gibbs     <- array(0,dim=c(p,p,nmc))
	E_gibbs     <- array(0,dim=c(p,p,nmc))
	delta_gibbs <- array(0,dim=c(p,1,nmc))


    v0_gibbs    <- rep(0,nmc)
	tau_update  <- array(0,dim=c(p,p,nmc)) 
    pi_update   <- rep(0,nmc)

    utriFlag <- upper.tri(matrix(NA,p,p)) ## Upper triangular flag


	for( gg in 1:(burnin+nmc)){##gg<-1

		########################
		#### Sample zhat_0  ####
		########################
		for( ii in 1:n){
			tmp  <- zsample(1,x[ii,],zhat[ii,],R_mc, delta_mc)
			zhat_mc[ii,tmp$ind0] <- as.vector( tmp$z0new ) 
			}


		########################
		#### Sample delta   ####
		########################

		zmax_mc <- rep(0,p) # variables-wise maximum among the censored observations
		zmin_mc <- rep(0,p) # variables-wise maximum among the censored observations
		for( jj in 1:p){
			zmax_mc[jj] <- max( zhat_mc[x[,jj]==0,jj] )
			zmin_mc[jj] <- min( zhat_mc[x[,jj]!=0,jj] )
		}

		if( sum(zmax_mc>zmin_mc)!=0  ){ warning("zmin has to be larger than zmax")}
		delta_mc <- runif(p,zmax_mc,zmin_mc)
		delta_mc[is.nan(delta_mc)] <- -Inf


		########################################
		#### Sample Sigma (Omega) and Graph ####
		########################################

        v1_mc = h * v0_mc   # (h in code) =  (h in paper )^2
        
        V0 = v0_mc * matrix(1, p, p); V1 = v1_mc * matrix(1, p, p)
        
        S    <- t(zhat_mc)%*%zhat_mc
        
		SSSL <- SSVS(S, n, Sig = R_mc, V0, V1, tau_mc, lambda, pijk_mc, burnin=0, nmc=1) 


		R_mc <- cov2cor( SSSL$Sig_save[,,1] )
		C_mc <- SSSL$C_save[,,1]
		E_mc <- SSSL$Z_save[,,1]
		tau_mc <- SSSL$T_save[,,1]	
		
		
        ########################################
        ####         Sample v0^2            ####
        ########################################
        
        
        pp <- p*(p-1)/2
        wijsq <- C_mc[utriFlag]^2
        
        v0rate  <- sum( wijsq/(h^(E_mc[utriFlag]) ) )/2 + v0beta
        v0shape <- pp/2 + v0alpha
        
        v0_mc   <- 1/rgamma(1, shape=v0shape, rate=v0rate)
        
        
                
        ########################################
        ####         Sample pijk_mc         ####
        ########################################
        
        ne1 <- sum(E_mc[utriFlag])
        pijk_mc <- rbeta(1,1+ne1,pp-ne1)
        
        
        
        
        
        
        
        

        
        
        
        
        
        
        
        
        
        
        
        
		if( gg > burnin){
			delta_gibbs[,,gg-burnin] <- delta_mc
			R_gibbs[,,gg-burnin] <- R_mc
			C_gibbs[,,gg-burnin] <- C_mc
			E_gibbs[,,gg-burnin] <- E_mc
            v0_gibbs[gg-burnin]  <- v0_mc
			tau_update[,,gg-burnin] <- tau_mc
            pi_update[gg-burnin] <- pijk_mc
            

        

            }

            if( TRUE ){
                if( gg%%1000 == 0 ){ print(gg);  print( Sys.time() )  }
            }

}


gibbsSample <- list(
			R_gibbs = R_gibbs,
			C_gibbs = C_gibbs,
			E_gibbs = E_gibbs,
            pi_gibbs = pi_update,
            v0_gibbs = v0_gibbs
			)


            
            
}










