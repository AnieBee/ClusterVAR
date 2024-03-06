
EMFunc <- function(Init,
                   IDNames,
                   Y,
                   X,
                   K,
                   N,
                   Tni_NPred,
                   qqq,
                   nDepVar,
                   NewPredictableObs,
                   LaggedPredictObs,
                   PredictableObsConc,
                   LaggedPredictObsConc,
                   PersEnd,
                   PersStart,
                   PersStartU_NPred,
                   PersEndU_NPred,
                   Covariates,
                   Conv,
                   it,
                   smallestClN,
                   SigmaIncrease)
{

  likelihood = rep(NA, it)
  iterationReset = rep(FALSE, it)
  # Have any likelihoods/posteriors been reset in this iteration (FALSE means no reset in this iteration)
  EMiteration = 0 # counts the iteration
  ratio = 1000
  lik = -10000000
  lowest.Likelihood = -9000

  FYZ = matrix(0, N, K) # prob(W_{i}| Sigma_k, A_K, B_K), Normal weighted by prior
  FZY = matrix(0, N, K) # posterior group membership (pi_{ik})
  A = Init$A
  B = Init$B
  Sigma = Init$Sigma # Sigma has to be checked for singularity, is checked for it at end of EMInit
  tau = Init$tau
  Lags = Init$Lags # contains switched Lags order to match lag order the clusters exhibited in EMInit
  U = array(NA, dim = c(nDepVar, PersEndU_NPred[[min(Lags)]][N], K))
  # Vector of u_{ikt}s # U is not of same length as Y, Y contains N many Lags*m pre-samples

  # Constraints on B #
  DimensionsBasedonConstraints = constraintsOnB(Covariates, K, N)
  Wk = array(0, dim = c(nDepVar, PersEnd[N], DimensionsBasedonConstraints$WkNumbVersions))
  # Wk is of dim(3) = K if B are unequal, if B are equal OR  individual-specific  dim(3) = 1

  while ( (EMiteration < it) & (ratio > Conv) )
  {#### E & M steps -------------

    EMiteration = EMiteration + 1

    ########  E-STEP -----------------------------------------

    ## Check: Sigma ----------------
    Wk = calculateW(Covariates = Covariates, K = K, Wk = Wk, Y = Y, B = B, X = X)

    U = calculateU(K = K, WkNumbVersions = DimensionsBasedonConstraints$WkNumbVersions,
                   PredictableObsConc = PredictableObsConc,
                   LaggedPredictObsConc = LaggedPredictObsConc, Tni_NPred = Tni_NPred,
                   U = U, Wk = Wk, A = A,
                   Lags = Lags, nDepVar = nDepVar)

    # calculate the FYZs -----------
    # all FYZs and FZYs are in log form to avoid underflow
    # E step, or dvmnorm could be implemented in Rcpp
    FYZ = calculateFYZ(K = K, N = N, FYZ = FYZ, U = U, 
                       PersStartU_NPred = PersStartU_NPred,
                       PersEndU_NPred = PersEndU_NPred,
                       nDepVar = nDepVar, Sigma = Sigma,
                       Lags = Lags)

    ## Check: likelihoods ----------
    FYZListCL = checkLikelihoodsNA(FYZ = FYZ, EMiteration = EMiteration)
    FYZ = FYZListCL$FYZ

    # calculate FZY
    calcPostList = calculatePosterior(N = N, K = K, tau = tau, FYZ = FYZ,
                                      lowest.Likelihood = lowest.Likelihood,
                                      FZY = FZY,
                                      UseFZY = TRUE)
    FZY = calcPostList$FZY
    # calculate FZY = posterior(pi_{ik})


    ### Checks --------------------
    ## Check for outliers
    FZYListCO = checkOutliers(K = K, FZY = FZY)

    ## Check posteriors are not NA ##
    FZYListCP = checkPosteriorsNA(FZY = FZYListCO$FZY, K = K)

    ## Check no component collapses onto a single point
    FZYListCCC = checkComponentsCollapsed(K = K, N = N, FZY = FZYListCP$FZY,
                                          smallestClN = smallestClN)
    FZY = FZYListCCC$FZY

    ########  M-STEP    ########
    ### calculate A -----------------
    A = calculateA(K = K, WkNumbVersions = DimensionsBasedonConstraints$WkNumbVersions,
                   N = N, Wk = Wk, NewPredictableObs = NewPredictableObs,
                   LaggedPredictObs = LaggedPredictObs, Lags = Lags, FZY = FZY,
                   A = A, nDepVar = nDepVar)
    # those places of A that will not be filled (because of lower lag number in some clusters)
    # will always be zero anyway because they are never filled, once A is calculated
    # in EMFunc, cannot be changed back

    ### calculate Sigma (S) ---------------------
    Sigma = calculateSigma(K = K, N = N, FZY = FZY, U = U, PersStartU_NPred = PersStartU_NPred,
                           PersEndU_NPred = PersEndU_NPred, Tni_NPred = Tni_NPred, Sigma = Sigma, Lags = Lags)
    SigmaList = checkSingularitySigma(nDepVar = nDepVar, K = K, Sigma = Sigma)
    Sigma = SigmaList$Sigma
    Sigma[ , , FZYListCCC$resetCl] = Sigma[ , , FZYListCCC$resetCl] + SigmaIncrease # Increase variance of components indicated by FZYListCCC

    
    ## Calculate B depending on Covariate constraint -------------
    B = calculateB(Covariates = Covariates, K = K, nDepVar = nDepVar, A = A,
                   Sigma = Sigma, N = N, NewPredictableObs = NewPredictableObs,
                   X = X, Y = Y, Lags = Lags, FZY = FZY,
                   qqq = qqq, B = B)

    ### Calculate tau -------------------
    tau = calculateTau(tau = tau, FZY = FZY, N = N)

    ### end of M-step ###

    #### Calculate current log likelihood (temp) and compare to old likelihood (lik) ---
    # calculating the denominator of pi _{ik} using updated tau
    calcPostList2 = calculatePosterior(N = N, K = K, tau = tau, FYZ = FYZ,
                                       lowest.Likelihood = lowest.Likelihood,
                                       FZY = FZY,
                                       UseFZY = FALSE)
    likelihood[EMiteration] = calcPostList2$logLikelihood

    ### Has EM been reset in this iteration? --------------
    iterationReset[EMiteration] =  FYZListCL$iterationReset |
      calcPostList$iterationReset |
      FZYListCP$iterationReset | FZYListCCC$iterationReset |
      SigmaList$iterationReset | FZYListCO$iterationReset |
      calcPostList2$iterationReset



    ## compare to old log likelihood (lik)
    ratio = calculateRatio(temp = likelihood[EMiteration], lik = lik,
                           EMiteration = EMiteration, Conv = Conv,
                           L2iterationReset = ifelse(EMiteration > 1,
                                                     iterationReset[EMiteration - 1] | iterationReset[EMiteration],
                                                     iterationReset[EMiteration]))
    lik = likelihood[EMiteration]
  } ### end of while loop (end of E & M steps)-------

  Converged = checkConvergence(Conv = Conv, ratio = ratio,
                               L2iterationReset = iterationReset[EMiteration - 1]
                               | iterationReset[EMiteration])
  nPara = calculateNPara(Lags = Lags, nDepVar = nDepVar, K = K,
                         BnumbVersions = DimensionsBasedonConstraints$BNumbVersions,
                         ncovariates = qqq)  
  Classification = apply(FZY, 1, which.max)
  last.lik = likelihood[EMiteration]
  
  ## claculate clTimepoints
  clTimepoints = as.vector(rep(0, K))
  for(j in 1:K)
  {
      for (i in 1:N)
      {
          clTimepoints[j] = clTimepoints[j] + (FZY[ i, j] * (Tni_NPred[[Lags[j]]][i]))
      }
  }##
  SC = calculateIC(ICType = "SC", Sigma = Sigma, Lags = Lags,
                   nDepVar = nDepVar, K = K, clTimepoints = clTimepoints, tau = tau)

  HQ = calculateIC(
      ICType = "HQ",
      Sigma = Sigma,
      Lags = Lags,
      nDepVar = nDepVar,
      K = K,
      clTimepoints = clTimepoints,
      tau = tau
  )
  
  BIC = calculateBIC(nPara = nPara, clTimepoints = clTimepoints, last.lik = last.lik)
  ICL = calculateICL(BIC = BIC, K = K, N = N, FZY = FZY, 
                     Classification = Classification)

  # Use ID names to return Classification and user knows what classification means
  Classification <- matrix(Classification, nrow = 1)
  colnames(Classification) <- IDNames

  colnames(B) = rownames(X) # Name every col in Array with Covariates Variable name
  rownames(B) = rownames(Y) # Name every row in B with Endogenous Variable Names
  colnames(A) = paste(rep(rownames(Y), max(Lags)), "_t-", sort(rep(1:max(Lags), length(rownames(Y)))), sep = "")
  rownames(A) = rownames(Y)

  #Intercept is the "reference" group in case categorical variables are included


  # ----- Return Output list -----
  outlist <- list(Converged = Converged,
                  VAR_coefficients = A,
                  Exogenous_coefficients = B,
                  Sigma = Sigma,
                  EMRepetitions = EMiteration,
                  last.loglik = last.lik,
                  nParameters = nPara,
                  LogLikelihood = likelihood[1:EMiteration],
                  Lags = Lags,
                  EMiterationReset = iterationReset[1:EMiteration],
                  Posterior_probabilities = t(FZY),
                  Classification = Classification,
                  Proportions = tau,
                  PredictableTimepoints = sum(clTimepoints),
                  SC = SC,
                  HQ = HQ,
                  BIC = BIC,
                  ICL = ICL)

  return(outlist)


} # end of EMfunc





