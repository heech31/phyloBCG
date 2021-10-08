synthData_from_ecdf_and_z <- function (comm, mar = 2, Sigma, n, seed = 10010, verbose = FALSE) 
{
    d <- ncol(comm)
    zratio <- apply(comm, MARGIN = mar, function(x) (sum(x == 
        0)/length(x)))
    maxabund <- apply(comm, MARGIN = mar, max)
    if (!is.null(seed)) {
        set.seed(seed)
    }
    normd <- mvtnorm::rmvnorm(n, mean = rep(0, d), sigma = Sigma)
    unif <- pnorm(normd)
    dat <- matrix(0, n, d)
    for (j in 1:d) {
        nzind <- which(unif[, j] > zratio[j])
        empf <- ecdf(comm[, j])
        ptm <- proc.time()
        for (k in 1:length(nzind)) {
            dat[nzind[k], j] <- qstepcdf(unif[nzind[k], j], empf, 
                interval = c(0, maxabund[j]))
        }
        if (verbose == TRUE) {
            cat("iteration = ", j, ": time = ", proc.time() - 
                ptm, "\n")
        }
    }
    return(list(dat=dat, z=normd) )
}




environment( synthData_from_ecdf_and_z ) <- asNamespace("SPRING")
assignInNamespace( "synthData_from_ecdf", synthData_from_ecdf_and_z, ns="SPRING" )





