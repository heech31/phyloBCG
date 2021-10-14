quantile.stepcdf <- function(prob, empf, interval, tol = 1e-3, maxiter = 100){
  # This function is from SPRING package
  # https://github.com/GraceYoon/SPRING
  ans <- c()
  # uniroot.all from rootSolve package was the fastest one.
  sol <- as.numeric(
    rootSolve::uniroot.all(function(x){empf(x)-prob}, interval = interval, tol = tol, maxiter = maxiter))
  
  if( length(sol) != 0 ){
    # in case quantile function is a step function, added the following step.
    if (prob <= empf(floor(sol))) {
      ans <- floor(sol)
    } else if (prob > empf(floor(sol))) {
      ans <- ceiling(sol)
    }
  }else{
    ans <- 0
  }
  return(ans)
}