calculateFYZ <- function(K, N, FYZ, U, PersStartU_NPred, PersEndU_NPred, nDepVar, Sigma, Lags)
{
    for(j in 1:K){  
        for(i in 1:N)
        {
            Uslice = t(U[ , (PersStartU_NPred[[Lags[j]]][i]):(PersEndU_NPred[[Lags[j]]][i]), j])
            FYZ[i, j] = sum(mvtnorm::dmvnorm(x = Uslice, 
                                    mean = as.matrix(rep(0, nDepVar), ncol = nDepVar),
                                    sigma = as.matrix(Sigma[, , j], ncol = nDepVar) , log = TRUE)) 
            # Sum of all log probabilities (from normal pdf) of columns of U, 
            # columns of x represent time points per individual 
        }
    } 
    
    invisible(FYZ)
}