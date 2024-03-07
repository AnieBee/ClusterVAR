calculateA <-
    function(K,
             WkNumbVersions,
             N,
             Wk,
             NewPredictableObs,
             LaggedPredictObs,
             Lags,
             FZY,
             A,
             nDepVar)
    {
        for (j in 1:K)
        {
            WkRunner = ifelse(WkNumbVersions == K, j, 1)
            
            AKnum = 0 # Sum of individual sums for A weighted by memb (tau), within j
            AKdenom = 0
            for (i in 1:N)
            {
                WK_Lagged = matrix(Wk[, LaggedPredictObs[[ Lags[j] ]][[i]], WkRunner], nrow = (nDepVar * Lags[j]))
                
                AKnum = AKnum + (FZY[i, j] * (Wk[, NewPredictableObs[[ Lags[j] ]][[i]], WkRunner] %*% t(WK_Lagged)))
                AKdenom = AKdenom + (FZY[i, j] * (WK_Lagged %*% t(WK_Lagged)))
            }
            A[, 1:(nDepVar * Lags[j]), j] = AKnum %*% MASS::ginv(AKdenom)
        }
        
        invisible(A)
        
    }