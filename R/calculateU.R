calculateU <- function(K, WkNumbVersions, PredictableObsConc, LaggedPredictObsConc, Tni_NPred, U, Wk, A, Lags, nDepVar)
{

    for(j in 1:K)
    {
        # Determine runner for W based on the constraints on B -------------------------
        WkRunner = ifelse(WkNumbVersions == K, j, 1) 

        # calculate U ------
        # If you would not index in A, you would not use all the Us you have at your disposal, even when lag number is smaller 
        U[ , 1:(sum(Tni_NPred[[ Lags[j] ]])), j] = matrix(Wk[ , PredictableObsConc[[ Lags[j] ]], WkRunner], nrow = nDepVar) - 
            ( A[ , 1:(nDepVar * Lags[j]), j] %*% matrix(Wk[ , LaggedPredictObsConc[[ Lags[j] ]], WkRunner], nrow = (nDepVar * Lags[j]) ))

    }
    
    invisible(U)
    
}