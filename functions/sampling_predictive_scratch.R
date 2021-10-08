# 
# is.symmetric.matrix( R_mc )
# is.symmetric.matrix( round(R_mc,9) )
# 
# R_mc <- ( R_mc + t(R_mc) )/2
# 
# is.positive.definite( R_mc )
# 
# x.new  <- x[1,]
# x <- x[-1,]
# 
# 
# 
# 
# eFx  <- apply(x,2,ecdf)
# eFxx <- Map(function(f,x) do.call(f, list(x)), eFx, alply(x,2)  )
# zhat <- matrix( unlist( lapply(eFxx,function(pr) qnorm( ( n/(n+1) )*pr ) ) ), n-1, p)
# maxabund <- apply(x, 2, max) # to restrict the search range of the solution (for uniroot.all)
# 
# # Flag indicating observed variables of x.new
# o.flag <- x.new!=0
# # Probability of Xj less than equal to xj
# F.new  <- ( n/(n+1) )*unlist( Map(function(f,x) do.call(f, list(x)), eFx, alply(x.new,1) ) )
# # Get z.o and z.t
# z.new  <- qnorm(F.new)
# z.new[!o.flag] <- 0
# 
# 
# z.new.loo <- rep(0,p)
# x.new.loo <- rep(0,p)
# 
# for( jj in 1:p){
#    
#   R.no.j       <- R_mc[-jj,-jj]
#   z.new.no.j   <- z.new[-jj]
#   o.flag.no.j  <- o.flag[-jj]
#   delta.t.no.j <- delta_mc[-jj][!o.flag.no.j]
#   # Partition correlation matrix
#   R.oo <- R.no.j[ o.flag.no.j, o.flag.no.j] # Variance-Covariance matrix of observed variables
#   R.ot <- R.no.j[ o.flag.no.j,!o.flag.no.j] # Covariance matrix between observed and truncated variables
#   R.tt <- R.no.j[!o.flag.no.j,!o.flag.no.j] # Variance-Covariance matrix of truncated variables
# 
#   iR.oo.ot <- solve(R.oo, R.ot) # invS11 * S10
#   c.mu    <- as.vector( crossprod( z.new.no.j[o.flag.no.j], iR.oo.ot) )
#   c.Sigma <- R.tt - crossprod(R.ot, iR.oo.ot)
#   c.Sigma <- (c.Sigma + t(c.Sigma))/2
#   
#   from.truncated.mvnorm <- rtmvnorm(1, mean=c.mu, sigma= c.Sigma, 
#            lower= rep(-Inf,length(c.mu)), upper = delta.t.no.j,
#            algorithm="gibbs",burn.in.sample=100)
#   
#   z.new.no.j[!o.flag.no.j] <- as.vector(from.truncated.mvnorm)
#   
#   mu.jj.given._jj <- sum( solve( R.no.j,R_mc[jj,-jj] )*z.new.no.j )
#   sd.jj.given._jj <- sqrt( R_mc[jj,jj] - sum( solve( R.no.j, R_mc[jj,-jj] ) * R_mc[jj,-jj] ) )
#   
#   z.new.loo[jj] <- rnorm(1,mu.jj.given._jj, sd.jj.given._jj)
#   z.new.loo[jj] <- ifelse( z.new.loo[jj]<delta_mc[jj], 0, z.new.loo[jj] )
#   
# }
# 
# plot( z.new.loo, z.new );abline(a=0,b=1,col=2)
# mean( (z.new.loo- z.new)^2 )
# #R_mc <- SigmaTrue
# #R_mc <- cor(zhat_mc)
# #R_mc <- solve( SigmaTrue )
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
# delta_mc[jj]
# mean( (z.new.loo-z.new)^2 )
# 
# z.new[jj]
# 
# 
# R_mc[,1]
# 
# 
# 
# 
# 
# for ( jj in 1:p ){
#   if( o.flag[[jj]] ){
#     # if the z.new_j is observed get the x.new_j using the estimated cdf
#     #nzind <- which(unif[, j] > zratio[j]) # to keep the zero ratio as the true data.
#     #empf <- ecdf(refdata[, j]) # empf is a cdf function. empf(c) = Pr(X <= c)
#     interval <- c(0, maxabund[jj])
#     x.new[jj] <- qstepcdf(u.new[jj], eFx[[jj]], interval = interval)
#   }
# }
# 
# 
# 
# 
# 
# crossprod(R.ot, iR.oo.ot)
# 
# t(R.ot)%*% solve(R.oo) %*% R.ot
# 
# cmu    <- as.vector( crossprod( z1, iS1101  ) )  # Conditional mean S01*invS11*z1
# cSigma <- S00 - S01%*%iS1101           # Conditional covariance S00 - S01*invS11*S10
# cSigma <- ( cSigma + t(cSigma) )/2
# lowerb <- rep(-Inf,p0) # lower limits for truncated normal
# upperb <- delta0 # upper limits for truncated normal
# 
# if(below){ # If truncated below
#   upper <- delta0
#   lower <- rep(-Inf,p0)
# }else{ # If truncated above
#   upper <- rep(Inf,p0)
#   lower <- delta0
# }
# 
# zt.new  <- tmvtnorm::rtmvnorm(1, mean=cmu, sigma=cSigma, upper=upper, lower=lower,
#                               algorithm = "gibbs", burn.in.samples=50, thinning= 1 )
# if( sum(is.na(zt.new))>0 ){ # If NA occurs
#   zt.new  <- tmvtnorm::rtmvnorm(1, mean=cmu, sigma=(cSigma+diag(0.1,p0)), upper=upper, lower=lower,
#                                 algorithm = "gibbs", burn.in.samples=50, thinning= 1 )
#   occur.fail <- 1
# }else{
#   occur.fail <- 0
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
# c( qnorm( c(eFx[[jj]]( x.new[jj] ), u.new[jj] )  ), z.new[jj] )
# 
# install.packages("rootSolve")
# 
# qstepcdf <- function(p, empf, interval, tol = 1e-3, maxiter = 100){
#   ans <- c()
#   # uniroot.all from rootSolve package was the fastest one.
#   sol <- as.numeric(rootSolve::uniroot.all(function(x){empf(x)-p}, interval = interval, tol = tol, maxiter = maxiter))
# 
#   # in case quantile function is a step function, added the following step.
#   if (p <= empf(floor(sol))) {
#     ans <- floor(sol)
#   } else if (p > empf(floor(sol))) {
#     ans <- ceiling(sol)
#   }
#   return(ans)
# }
# 
# 
