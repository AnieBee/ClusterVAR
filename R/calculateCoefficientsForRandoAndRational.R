

calculateCoefficientsForRandoAndRational <- function(Covariates,
                                                     K,
                                                     N,
                                                     nDepVar,
                                                     qqq,
                                                     Lag,
                                                     PersEnd,
                                                     PersStart,
                                                     Y,
                                                     X,
                                                     NewPredictableObs,
                                                     LaggedPredictObs) {


  DimensionsBasedonConstraints = constraintsOnB(Covariates, K, N)
  BIndividual = array(NA, dim = c(nDepVar, qqq, N))
  CoefficientsForRandoAndRational = array(NA, dim = c((nDepVar * nDepVar * Lag) +
                                                        (DimensionsBasedonConstraints$ClusterOnB * (nDepVar * qqq)),
                                                      N))
  WIndividual = array(NA, dim = c(nDepVar, PersEnd[N]))
  for(i in 1:N)
  {
     XPerson = X[ , PersStart[i]:PersEnd[i], drop = FALSE]
    BIndividual[ , , i] = (Y[ , PersStart[i]:PersEnd[i]] %*% t(XPerson)) %*%
      ginv(XPerson %*% t(XPerson))
    WIndividual[ , (PersStart[i]):(PersEnd[i])] = Y[ , (PersStart[i]):(PersEnd[i])] - (BIndividual[ , , i] %*% XPerson)
    
      LaggedWK = matrix(WIndividual[ , LaggedPredictObs[[i]], drop = FALSE], nrow = (nDepVar * Lag))
      AKn = WIndividual[ , NewPredictableObs[[i]], drop = FALSE] %*% t(LaggedWK)
      AKd = LaggedWK %*% t(LaggedWK)
    CoefficientsForRandoAndRational[1:(nDepVar * nDepVar * Lag), i] = as.vector(AKn %*% ginv(AKd)) # AIndividual # gets read in left to right
  }
  if(DimensionsBasedonConstraints$ClusterOnB)
  {# Add B to these coefficients on which the initial clustering solutions are built
    CoefficientsForRandoAndRational[(1 + (nDepVar * nDepVar * Lag)):((nDepVar * qqq) + (nDepVar * nDepVar * Lag)),  ] = as.vector(BIndividual)
    # gets read in top to bottom
  }

  return(CoefficientsForRandoAndRational)

} # EoF
