calculateU <- function(K, WkNumbVersions, N, NewPredictableObs, LaggedPredictObs,  Tni_NPred, U, Wk, A, Lags, nDepVar)
{

    for(j in 1:K)
    {
        # Determine runner for W based on the constraints on B -------------------------
        WkRunner = ifelse(WkNumbVersions == K, j, 1) 
        PredictableObsConc = do.call(c, NewPredictableObs[[ Lags[j] ]]) 
        LaggedPredictObsConc = do.call(c, LaggedPredictObs[[ Lags[j] ]]) 
        
        # calculate U ------
        # If you would not index in A, you would not use all the Us you have at your disposal, even when lag number is smaller 
        U[ , 1:(sum(Tni_NPred[[ Lags[j] ]])), j] = matrix(Wk[ , PredictableObsConc, WkRunner], nrow = nDepVar) - 
            ( A[ , 1:(nDepVar * Lags[j]), j] %*% matrix(Wk[ , LaggedPredictObsConc, WkRunner], nrow = (nDepVar * Lags[j]) ))

    }
    
    invisible(U)
    
}