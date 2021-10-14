SSVS  = function(S, n, Sig, V0, V1, tau, lambda, pijk) {

  ## Scale-up Stochastic Structural Learning (SSSL, Wang 2014) to a concentration graph models
  ## **NOTE** This function implement ONE Gibbs update
  ## Input:
  ##   S: p x p gram matrix of previous step:  X'*X if X is a n x p matrix of random sample
  ##   n: sample size
  ##   Sig: Sigma correlation matrix from the previous step
  ##   V0:  p x p matrix of the small variance  components;  V0(j,k) = v0 for all j,k , v0 is from the previous step
  ##   V1:  p x p matrix of the large variance  components;  V1(j,k) = h*v0 for all j,k
  ##   tau: p x p matrix of the variance components (spike and slab) corresponding to C(j,k) from the previous stp
  ##   lambda: Hyperparameter for the diagonal element. Usually set to be 1;
  ##   pijk: p x p matrix of the edge inclusion probabilities from the previous step

  ## Update order: 
  
  ## Output:
  ##  Sig_new: p x p  saved covariance matrices
  ##  C_new: p x p saved concentration matrices
  ##  Z_new: p x p saved graph indicator matrices
  ##  T_new: p x p saved variance component matrix (tau)
  
  
  p = nrow(S)
  
  Z =  matrix(0, p, p) # Adjacency matrix 
  C =  matrix(0, p, p) # Concentration matrix
  
  ind_noi_all = matrix(0, p - 1, p)
  
  for (i in 1:p) { #(p -1) x p index matrix without the ith index
	ind_noi_all[, i] = setdiff(1:p,i)
  }



    ### sample Sig and C = inv(Sig)
    for (i in 1:p) {
      ind_noi = ind_noi_all[, i]
      tau_temp = tau[ind_noi, i]
      Sig11 = Sig[ind_noi, ind_noi]
      Sig12 = Sig[ind_noi, i]
      
      
      invC11 = Sig11 - Sig12 %*% t(Sig12) / Sig[i, i]
      Ci = (S[i, i] + lambda) * invC11 + diag(1. / tau_temp)
      Ci = (Ci + t(Ci)) / 2
      Ci_chol = chol(Ci)
      
      mu_i = -solve( Ci_chol,  solve(t(Ci_chol), S[ind_noi, i]) )
      beta = mu_i + solve(Ci_chol, rnorm(p - 1))
      
      C[ind_noi, i] = beta
      C[i, ind_noi] = beta
      
      
      a_gam = 0.5 * n + 1
      b_gam = (S[i, i] + lambda) * 0.5
      gam = rgamma(1, shape = a_gam, rate = b_gam)
      c = t(beta) %*% invC11 %*% beta
      C[i, i] = gam + c
      
      ## Below updating Covariance matrix according to one-column change of precision matrix
      invC11beta = invC11 %*% beta
      Sig[ind_noi, ind_noi] = invC11 + invC11beta %*% t(invC11beta) / gam
      Sig12 = -invC11beta / gam
      Sig[ind_noi, i] = Sig12
      Sig[i, ind_noi] = Sig12
      Sig[i, i] = 1 / gam
      
      if( length(pijk)>1 ){ # pijk from the previous step
	      pii     = pijk[ind_noi,i] # 
  	    } else{
	  	  pii     = pijk # constant pijk (Oracle case)
	      }
      
      
      v0 = V0[ind_noi, i]
      v1 = V1[ind_noi, i]
      
      w1 = -0.5 * log(v0) - 0.5 * beta ^ 2 / v0 + log(1 - pii)
      w2 = -0.5 * log(v1) - 0.5 * beta ^ 2 / v1 + log(pii)
      
      w_max = pmax(w1, w2)
      w = exp(w2 - w_max) / rowSums(exp(cbind(w1 - w_max, w2 - w_max)))
      
      z = (runif(p - 1) < w)
      v = v0
      v[z] = v1[z]
      
      
      tau[ind_noi, i] = v # Update variance components (spike or slab)
      tau[i, ind_noi] = v # Update variance components (spike or slab)
      
      Z[ind_noi, i] = z # Update edge indicators
      Z[i, ind_noi] = z # Update edge indicators
      
  }
  
  return (list( Sig_new = Sig, C_new = C, Z_new = Z, T_new = tau)  )

  }


