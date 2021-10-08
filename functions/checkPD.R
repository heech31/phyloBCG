library(matrixcalc) #is.positive.definite
checkSymmetricPositiveDefiniteCustom <- function (x, name = "sigma")
{
    if (!isSymmetric(x, tol = sqrt(.Machine$double.eps))) {
        stop(sprintf("%s must be a symmetric matrix", name))
    }
    if (NROW(x) != NCOL(x)) {
        stop(sprintf("%s must be a square matrix", name))
    }
    if (any(diag(x) <= 0)) {
        stop(sprintf("%s all diagonal elements must be positive", 
            name))
    }
    if ( !is.positive.definite(x, tol=1e-8) ) {
        stop(sprintf("%s must be positive definite", name))
    }
}
