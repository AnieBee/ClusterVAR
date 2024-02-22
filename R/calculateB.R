calculateB <- function(Covariates, K, nDepVar, A, Sigma, N, NewPredictableObs, X, Y, Lags, FZY, qqq, B)
{
    if(Covariates == "equal-within-clusters"){# B(0): equal-within-clusters
        for(j in 1:K)
        {
            Bnum = 0
            Bdenom = 0
            Delta = cbind(diag(nDepVar), (-1) * A[ , 1:(nDepVar * Lags[j]), j])
            NewCovaMat = t(Delta) %*% ginv(Sigma[, , j]) %*% Delta
            for(i in 1:N)
            {
                Bn = 0
                Bd = 0
                for(trunner in c(NewPredictableObs[[ Lags[j] ]][[i]]) )
                { 
                    XtildaKron <- kronecker(t(X[ , (trunner):(trunner - Lags[j]), drop = FALSE]),
                                            diag(nDepVar)) # drop = FALSE is needed here in case of q = 1
                    Bd = Bd + ( t(XtildaKron) %*% NewCovaMat %*% XtildaKron ) 
                    Bn = Bn + ( t(XtildaKron) %*% NewCovaMat %*% 
                                     as.vector(Y[ , (trunner):(trunner - Lags[j]), drop = FALSE]) )
                }
                Bnum = Bnum + (FZY[ i, j] * Bn)
                Bdenom = Bdenom + (FZY[ i, j] * Bd)
            }
            BasVec = ginv(Bdenom) %*% Bnum
            B[, , j] = matrix(BasVec, nrow = nDepVar, ncol = qqq, byrow = FALSE)
        }
    }
    
    invisible(B)
    
}

