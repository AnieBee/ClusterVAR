calculateSigma <- function(K, N, FZY, U, PersStartU_NPred, PersEndU_NPred, Tni_NPred, Sigma, Lags)
{
    for(j in 1:K)
    {
        Snum = 0
        Sdenom = 0
        for(i in 1:N)
        { # there is a runner inside U for every individual i so the U multiplication result can be weighted by the posterior pi (memb)
            Uslice = U[ , (PersStartU_NPred[[Lags[j]]][i]):(PersEndU_NPred[[Lags[j]]][i]), j]
            Snum = Snum + (FZY[ i, j] * ( Uslice %*% t(Uslice) ) )
            Sdenom = Sdenom + (FZY[ i, j] * (Tni_NPred[[Lags[j]]][i]) )
        }
        Sigma[, , j] = Snum / Sdenom # Sdenom is integer
    }
    
    invisible(Sigma)
    
}