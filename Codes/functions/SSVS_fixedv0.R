SSVS  = function(S, n, Sig, V0, V1, tau, lambda, piij, burnin, nmc) {

  ## Almost verbatim translation from BayesGGM_SSVS_FixedV0V1.m
  ## BayesGGM_SSVS_FixedV0V1 fits SSSL (Wang 2014) to a concentration graph models
  
  ## Input:
  ##   S: p x p cross product matrix:  Y*Y' if Y is a p x n matrix of random sample
  ##   n: sample size
  ##   Sig: initial guess of Sigma
  ##   V0:  p x p matrix of the small variance  components;  V0(i,j) is the small variance of the precision element C(i,j)
  ##   V1:  p x p matrix of the large variance  components;  V0(i,j) is the small variance of the precision element C(i,j)
  ##   tau: p x p matrix of the variance components corresponding to C(i,j)
  ##   lambda: hyperparameter for the diagonal element Usually set to be 1;
  ##   piij: p x p matrix of the prior marginal edge "inclusion probabilities"
  ##   burnin: # of discarded samples
  ##   nmc: # of saved samples
  
  ## Output:
  ##  Sig_save: p x p x nmc  saved covariance matrices
  ##  C_save: p x p x nmc saved concentration matrices
  ##  Z_save: p x p x nmc saved graph indicator matrices
  ##  T_save: p x p x nmc saved variance component matrix (tau)
  
  C = solve(Sig)
  
  p = nrow(S)
  

  C_save   = array(0, c(p, p, nmc))
  Sig_save = C_save
  Z_save   = C_save
  T_save   = C_save
    
  Z =  matrix(1, p, p)
  
  
  ind_noi_all = matrix(0, p - 1, p)
  
  for (i in 1:p) {
	ind_noi_all[, i] = setdiff(1:p,i)
  }

  
  
  pii_RB = matrix(0, p, p)
  
  pii_mat = matrix(0, p, p)
  
  
  for (iter in 1:(burnin + nmc)) {
    if (iter %% 2000 == 0) {
      print(paste("iter = ", iter, " nedge = ", (sum(Z) - p) / 2, sep = ""))
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
      
      mu_i = -solve(Ci_chol, solve(t(Ci_chol), S[ind_noi, i]))
      
      beta = mu_i + solve(Ci_chol, rnorm(p - 1))
      
      
      
      C[ind_noi, i] = beta
      
      C[i, ind_noi] = beta
      
      
      a_gam = 0.5 * n + 1
      
      b_gam = (S[i, i] + lambda) * 0.5
      
      gam = rgamma(1, a_gam, b_gam)
      
      
      c = t(beta) %*% invC11 %*% beta
      
      C[i, i] = gam + c
      
      
      
      
      
      
      ## Below updating Covariance matrix according to one-column change of precision matrix
      invC11beta = invC11 %*% beta
      
      
      Sig[ind_noi, ind_noi] = invC11 + invC11beta %*% t(invC11beta) / gam
      
      Sig12 = -invC11beta / gam
      
      Sig[ind_noi, i] = Sig12
      
      Sig[i, ind_noi] = Sig12
      
      Sig[i, i] = 1 / gam
      


      
      if( length(piij)>1 ){
	      pii     = piij[ind_noi,i]
	  } else{
	  	  pii     = piij
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
      
      
      
      
      
      pii_mat[ind_noi, i] = w
      
      
      tau[ind_noi, i] = v
      
      tau[i, ind_noi] = v
      
      
      Z[ind_noi, i] = z
      
      Z[i, ind_noi] = z
      
      
      
      
    }
    
    
    
    
    
    if (iter > burnin) {
      pii_RB = pii_RB + pii_mat / nmc
      
      Sig_save[, , iter - burnin] = Sig
      
      C_save[, , iter - burnin] = C
      
      Z_save[, , iter - burnin] = Z

      T_save[, , iter - burnin] = tau
    }
    
    
    
  }
  
  return (list(
    Sig_save = Sig_save,
    C_save = C_save,
    Z_save = Z_save,
    T_save = T_save    )
    )
}


