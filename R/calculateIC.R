calculateIC <- function(ICType, Sigma, Lags, nDepVar, K, N, FZY, Tni, tau)
    # clusterTimePoints = T - k but extended to the cluster case
    # Calculate time-series IC for all clusters for a certain lag number (YWResidualCovariance, clusterTimePoints and LagsForP)
    # These are not global information criteria, they can only be used to compare timeseries models across a fixed number of clusters,
    # not across a different number of clusters
{

    ParasCl = as.vector(rep(0, K))
    clTimepoints = as.vector(rep(0, K))
    penaltyTerm = as.vector(rep(0, K))
    clIC = rep(0, K)
    
    for(j in 1:K)
    {
        # ParasCl[j] = calculateNPara(Lags = Lags[j], nDepVar = nDepVar, K = 1,
        #                             BnumbVersions = ifelse(BnumbVersions == 1,
        #                                                           tau[j],
        #                                                           1),
        #                             ncovariates = ncovariates)  + ((K - 1) / K) # includes all paras except tau, so + ((K - 1) / K)
        ParasCl[j] = Lags[j] * (nDepVar * nDepVar) # only include time-series parameters
        for (i in 1:N)
        {
            clTimepoints[j] = clTimepoints[j] + (FZY[ i, j] * (Tni[[Lags[j]]][i]))
        }
        
        penaltyTerm[j] = switch(ICType, 
                                 "HQ" = log(log(clTimepoints[j])), # is called HQ(n) in VARselect from vars package
                                 "SC" = log(clTimepoints[j]),
                                 #"AIC" = 1
                                )
        # AIC not offered anymore to users, they would get confused why it is so different to the global BIC, 
        # this is the AIC based on the determinant of sigma (the time-series AIC), the global BIC is based on the log likelihood
        
        ### Use ParasCl, clTimepoints and penaltyTerm to calculate clIC ###
        clIC[j] = tau[j] * ( log(det(Sigma[ , , j])) + ((2 * ParasCl[j] * penaltyTerm[j]) / clTimepoints[j]) )
        
    }
   
    
    invisible(sum(clIC))
    
}