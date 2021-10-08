# # Pivotal discrepancy measure (Biometrics, Yuan and Johnson, 2012)
# # for PhyloBCG
# # For truncated normal distributions, we use marginal cdfs as 
# # a pivotal quantity (Tibshirani et al., 2015)
# # Exact Post-Selection Inference for Sequential Regression Procedures (JASA)
# 
# PDM <- function(x, z, R, delta, nK = 4, nL = 11 ){
#   # data : needed for the locations of truncated variables
#   # z    : The sth Gibbs sample of Gaussian data
#   #        observed variables are the same across s=1,...,S
#   #        truncated variables are the sth Gibbs sample from truncated Gaussian
#   # R    : The sth Gibbs sample of Guassian correlation matrix
#   
#   n <- nrow(z) # Sample size
#   p <- ncol(z) # Number of variables
#   # list (n) that will store pivotal quantities and fitted value of zs
#   pdm    <- rep(NA,p) # PDM of p-variables
#   # Pivotal quantities and 
#   pivots <- fitted.z.obs <- matrix(NA,n,p)
# 
#   for( ii in 1:n){
#     i.obs.flag <- x[ii,] != 0 # flag of observed variables of ith obs
#     z.obs <- z[ii,  i.obs.flag] # observed Gaussian variables z.o
#     z.trc <- z[ii, !i.obs.flag] # truncated Gaussian variables z.t
#     lower <- delta[i.obs.flag] # Lower limit of observed Gaussians
#     
#     R.to  <- R[!i.obs.flag,  i.obs.flag] # Cov(z.t,z.o)
#     R.tt  <- R[!i.obs.flag, !i.obs.flag] # Var(z.t)
#     R.oo  <- R[ i.obs.flag,  i.obs.flag] # Var(z.o)
#     iR.tt.to <- solve(R.tt,R.to)         # {Var(z.t)}^{-1} Cov(z.t,z.o)
#     
#     # Conditioning on truncated z.t
#     mu.io    <- as.vector( crossprod(z.trc, iR.tt.to ) ) # Conditional mean
#     Delta.io <- R.oo - t(R.to)%*%iR.tt.to # Conditional variance
#     
#     if( !is.symmetric.matrix(Delta.io)){# If it is not symm, because of numerical error
#       Delta.io <- (Delta.io + t(Delta.io))/2 # Make it symmetric
#     }
#     
#     # Marginal standard deviations
#     sigmas <- sqrt( diag(Delta.io) )
#     # Denominator (Normalizing constant) of left-truncated normal cdf
#     denom  <- 1 - pnorm( (lower - mu.io)/sigmas )
#     # Numerator of left-truncated normal cdf
#     num    <- pnorm( (z.obs - mu.io)/sigmas ) - pnorm( (lower - mu.io)/sigmas )
#     
#     # Left-truncated normal cdf is a pivotal quantity
#     pivot  <- num/denom
#     pivots[ii, i.obs.flag] <- pivot
#     fitted.z.obs[ii, i.obs.flag] <- mu.io
#   }
#   # Group pivots in terms of fitted values (mu.io)
#   # and within each group partition, assign pivots
#   # to probability cells. See Section 2.3 of Yuan and Johnson (2012)
#   cutt.prop.fitted <- seq(0,1,length=2+nK-1)[-c(1,2+nK-1)]
#   cutt.prop.pivot  <- seq(0,1,length=2+nL-1)[-c(1,2+nL-1)]
#   
#   for( jj in 1:p){
#     # fitted values of jth variable and take only observed subjects
#     fitted.jth <- na.exclude(fitted.z.obs[,jj])
#     # Get quantiles of fitted.jth
#     cutt.point.fitted <- quantile( fitted.jth , cutt.prop.fitted)
#     # Since pivots follow Unif(0,1), 
#     # cutt points are the same as the cutt.prop.pivot
#     cutt.point.pivot  <- cutt.prop.pivot
# 
#     # Assign group labels
#     Kpartition.label <- rep(nK,length(fitted.jth))
#     for(kk in 1:(nK-1)){
#       Kpartition.label <- Kpartition.label - 1*(fitted.jth<cutt.point.fitted[kk])
#       }#End of kk
#     # Number of pivots in each group
#     nks <- as.vector( table( Kpartition.label ) )
#   
#     # discrepancy measure for each group
#     dks <- rep(0,nK)
#   
#     for( kk in 1:nK){
#       # pivots of jth variable and take only observed ones
#       pivot.jth <- na.exclude(pivots[,jj])[Kpartition.label==kk]
#       Lpartition.label <- rep(nL,nks[kk])
#       # Get cell label
#       for(ll in 1:(nL-1)){
#         Lpartition.label <- Lpartition.label - 1*(pivot.jth<cutt.point.pivot[ll])
#         }#End of ll
#     
#       # Observed counts of each cell
#       Ol <- table( factor(Lpartition.label, ordered=TRUE, levels = 1:nL) )
#       # Probabilities of each cell
#       pl <- diff( c(0,cutt.prop.pivot,1) )
#     
#       # Discrepancy measures for each group 
#       dks[kk] <- sum(  ( (Ol - nks[kk]*pl)/sqrt( nks[kk]*pl ) )^2  )
#       }# End of kk
#   
#     pdm[jj] <- sum( dks ) # Pivotal discrepancy measure of the jth variable
# 
#   }#End of jj
#   return(pdm)
#   
# }#End of function
# 
# 
# # From Yuan and Johnson (2012):
# # P(PDM_d > t) <= min( 1, potential.p.value)
# # Johnson (2007) states that
# # For example, if t represents the 0.99 quantile from G, 
# # it follows for large n that the probability that 
# # the S(.8n) > t is less than or equal to 0.05. 
# # In other words, a finding that the .8 quantile from 
# # the posterior distribution of pivotal quantities exceeds 
# # the .99 quantile from the nominal distribution implies that 
# # the PPP p value is less than .05.
# 
# 
# p.values <- function(PDM, nK=4, nL=11, ref.quantile=0.95){
# 
#   logical.for.p <- TRUE
#   r <- nmc-10  # Order index
#   t <- qchisq(ref.quantile,df=nK*(nL-1)) # Ref distribution quantile
#   ubound <- 1 # Initialize the upper bound of p.value
#   while(logical.for.p){
#     ubound.old <- ubound
#     # Upper-bound of p-value
#     ubound <- nmc*pchisq( t, df=nK*(nL-1), lower.tail = FALSE)/(nmc-r+1)
#     ubound <- min( c(1, ubound) )
#     # It this condition is true then the p.value is less than upper-bound
#     logical.for.p <- (quantile(PDM,r/nmc) > t)
#     
#     r <- r - 1
#   }
#   # While loop will stop when the condition is not met
#   # (this implies that the p.value is not smaller than 
#   # the current upper-bound (ubound)
#   # so we set the p.value = (ubound.old + ubound)/2
#   p.val <- (ubound.old + ubound)/2
#   return( p.val )
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Scratches
# 
# # apply(PDMs, 1, p.values, nK=3, nL=11, ref.quantile=0.95)
# # PDMs <- matrix(NA,p,nmc)
# # for( ss in 1:nmc){
# #   PDMs[,ss] <- PDM(x,gibbsSample$z_gibbs[,,ss],gibbsSample$R_gibbs[,,ss])
# # }
# # 
# # 
# # hist( PDMs[1,] )
# # 
# # qchisq(0.95,df=nK*(nL-1))
# # hist( PDM,prob=TRUE,30)
# # lines(seq(300,800,l=1000), dchisq( seq(300,800,l=1000), df=nK*(nL-1)) )
# # 
# # abline( v= qchisq(0.99,df=nK*(nL-1)), col=2)
# # mean( PDM>qchisq(0.99,df=nK*(nL-1)) )
# # 
# # 
# # r <- 4596
# # 
# # t <- qchisq(0.99,df=nK*(nL-1))
# # 
# # # log.terms <- 
# # #   lchoose(nmc, r:nmc) +
# # #   (r:nmc)*pchisq( t, df=nK*(nL-1), log.p = TRUE) +
# # #   (nmc-r:nmc)*pchisq( t, df=nK*(nL-1), log.p = TRUE, lower.tail = FALSE)
# # # 1 - sum( exp(log.terms) )
# # 
# # nmc*pchisq( t, df=nK*(nL-1), lower.tail = FALSE)/(nmc-r+1)
# # 
# # mean( PDM>t ) > (1-r/nmc)
# 
#   
#   
#   
# #}
# 
