

LCVARclust <- function(Data,
                       yVars,
                       Time,
                       ID,
                       xContinuous = NULL,
                       xFactor = NULL,
                       Covariates = "equal-within-clusters", # "equal-accros-clusters", "individual-specific"), #Anja: so far only equal-within-clusters is implemented
                       Clusters,
                       Lags,
                       smallestClN = 3,
                       RndSeed = 3,
                       Rand = 1,
                       Rational = TRUE,
                       Initialization = NULL,
                       SigmaIncrease = 10,
                       it = 25,
                       Conv = 1e-06,
                       pbar = TRUE,
                       ...)

# Data = data frame containing one or several columns for: yVars (position of columns containing endogenous VAR process variables in dataframe Data),
#   time point (position of column is indicated as integer with Time),
#   ID (integer giving position of column containing ID, a different ID variable for every participant),
#   continious exogenous variables (pos. indicated by xContinuous) and categorical exogenous variables (position of column is given by xFactor),
#   Initialization (gives position of a column that contains a guess at participants clustermembership)
# Rand = integer, how many random initialization should be carried out
# it = number of maximum iterations of EM algorithm
# Conv = Convergence criterion, accompanying publication


{

  # ------ Computing Some Aux Vars ------
  LowestLag = min(Lags)
  HighestLag = max(Lags)


  # ------ Input Checks ------

  # Cluster Search Sequence
  if(missing(Clusters)) stop("Specify the sequence of number of clusters to search.")

  # Lag Search Sequence
  if(missing(Lags)) stop("Specify the sequence of number of lags to search.")
  if(!all(Lags == round(Lags))) stop("Lags need to be specified as integers.")
  if(!all((Lags[-1] - Lags[-length(Lags)]==1))) stop("Lags need to be specified as subsequent integers.")

  # Set random seed, if specified
  if(!missing(RndSeed)) set.seed(RndSeed)

  # Checks
  stopifnot( ! (duplicated(c(ID, xFactor, xContinuous, Initialization))) ) # Evaluate there is no overlap
  stopifnot(LowestLag & HighestLag)
  stopifnot(HighestLag >= LowestLag)
  stopifnot(smallestClN > 1) # Smallest clustersize that is allowed
  # smallestClN is used in checkComponentsCollapsed
  stopifnot(Clusters > 0)
  stopifnot(length(yVars) > 1) # So far only multivariate time-series are implemented
  stopifnot(all(is.numeric(Clusters)))
  stopifnot(!any(duplicated(Clusters)))
  Clusters = Clusters[order(Clusters)]


  call <- match.call()

  # FZYcriterion = 1e-8
  PreviousSol = TRUE # Use solution of previous Lags as a start

  # Y has to be numeric, X can be numeric or factor
  # Y is data of form m \times (sum^N (nObs) )
  # ID = has to be factor
  ##### Preprocessing of Data Set #####--------------------
  Data = as.data.frame(Data)
  Data = Data[order(Data[ , ID], Data[ , Time]), ]
  # order Data according to ID, make sure an individual's observations occur one after another with first obs first, second second etc
  # observations have to occur ascending in time

  
  ### NEW------------ Insert rowwise deletion of NAs here------------------
  
  
  # Endogenous Variables #-------------------
  Y = t(as.matrix(Data[ , yVars]))
  nDepVar = dim(Y)[1]

  # Exogenous Variables #-----------------------------
  X = createX(YLength = dim(Y)[2], xFactor = xFactor, Data = Data, xContinuous = xContinuous)
  qqq = dim(X)[1] # Number of covariates variables, including intercept (q)


  # ID indicator Variable & TS length for every person #------------------------------
  Pers = as.factor(Data[ , ID])
  pers = unique(Pers)
  N = length(unique(Pers)) # length(levels(Pers))

  nObs = rep(0, N)
  for (i in 1:N)
  {# determine the total number of Observations (nObs) per pers
    nObs[i] = length(which(Pers == pers[i]))
  }

  stopifnot(identical(rep(pers, c(nObs)), as.factor(Data[ , ID]))) # should be superflous is a check nObs is correct
  stopifnot(nObs > (10 + HighestLag)) # check all indivudals have 11 more obs than lags
  stopifnot(is.numeric(Y))
  stopifnot(N >= (smallestClN * max(Clusters)))# check that for the highest number of clusters requested
  # there can be at least the smallest permissable number of people per cluster

  # Whole Sample of Obs
  PersStart = cumsum(c(0, nObs[-length(nObs)])) + 1 # Start of individual time series for every pers in Y or W
  PersEnd = cumsum(nObs) # end of individual time series of every pers, last obs of every pers in Y or W

  ### Runners that depend on the number of Lags ### ----
  PersPDiffStart = vector("list", HighestLag) # from 1:HighestLag but only LowestLag:HighestLag elements are filled
  Tni = vector("list", HighestLag)
  PersStartU = vector("list", HighestLag)
  PersEndU = vector("list", HighestLag)

  for (lagRunner in LowestLag:HighestLag)
  {
    PersPDiffStart[[lagRunner]] = PersStart + lagRunner # start of individual time series when first P observations
    # of every person are taken away, where p is the lag order of the VAR(p) process
    Tni[[lagRunner]] =  nObs - lagRunner  # length of U, lenght of sum of time series with presample of
    # first P observations removed, length of obs per pers in U
    # With presample of first Lags-obs removed for every individual (needed for U)
    PersStartU[[lagRunner]] = cumsum(c(0, nObs[-length(nObs)] - lagRunner)) + 1 # Start of individual time series for every pers in
    # U ( = presample is removed), first obs of a pers in U
    PersEndU[[lagRunner]] = cumsum(nObs - lagRunner) # Last obs of a pers in U
  }

  ### NEW------------ Insert determination of whether an observation in Y can be predicted (PredictableObs) ------------------
  PredictableObs = vector("list", HighestLag)  # from 1:HighestLag but only LowestLag:HighestLag elements are filled
  for (lagRunner in LowestLag:HighestLag) 
  {
      PredictableObs[[lagRunner]] = rbind( rep(1, dim(Y)[2]) , Pers)
      # assign zero to first column for non-predictable values
      PredictableObs[[lagRunner]][1, c(22, 32, 42)] = 0
      PredictableObs[[lagRunner]][1, c(PersStart)] = 0
      if(lagRunner > 1) PredictableObs[[lagRunner]][1, c(PersStart + 1)] = 0
      if(lagRunner > 2) PredictableObs[[lagRunner]][1, c(PersStart + 2)] = 0
      
  }
 ###### End of determining PredictableObs
  
  Tni_NPred = vector("list", HighestLag)  # from 1:HighestLag but only LowestLag:HighestLag elements are filled
  PersStartU_NPred = vector("list", HighestLag)
  PersEndU_NPred = vector("list", HighestLag)
  for (lagRunner in LowestLag:HighestLag) 
  {
      Tni_NPred[[lagRunner]] = rep(0, N)
      nPredObsPerPerson = PredictableObs[[lagRunner]][2, PredictableObs[[lagRunner]][1, ] == 1]
      for (i in 1:N)
      {# determine the total number of predictable Observations per person (Tni_NPred) differs across lag numbers
          Tni_NPred[[lagRunner]][i] = length(which(nPredObsPerPerson == pers[i]))
          ##----- NEW: do lags still have to be removed from Tni, basically this should be taken care of in PredictableObs already
          
          # Tni_NPred[[lagRunner]] # length of U, lenght of sum of time series with presample of
          # first P observations removed, length of obs per pers in U
          # With presample of first Lags-obs removed for every individual (needed for U)
      }
      ## Change PredictableObs so that it contains the position of the observations in Y, and W, and X that can be predicted
      PredictableObs[[lagRunner]] = which(PredictableObs[[lagRunner]][1,] ==  1)
     
      PersStartU_NPred[[lagRunner]] = cumsum(c(0, Tni_NPred[[lagRunner]][-length(Tni_NPred[[lagRunner]])])) + 1 # Start of individual time series for every pers in
      # U ( = presample is removed), first obs of a pers in U
      PersEndU_NPred[[lagRunner]] = cumsum(Tni_NPred[[lagRunner]]) # Last obs of a pers in U
      
  }
  
  ##----- NEW: PersStartU_NPred, PersEndU_NPred, Tni_NPred are new, are the old versions still needed anywhere? Also, PersPDiffStart is not 
  # needed anymore I think, because PredictableObs already removes all the first lag-many obs of every person so PersStart in combination with PredictableObs can be used 
  ##----- NEW: Add some kind of check that each person has enough predictable observations to estimate the model
  
 
  ### End ------------
 
  # Create val.init, a list containing memb, a vector ordered by ID that contains the cluster membership initialization ----
  if ( ! is.null(Initialization))
  {
    val.init = list(memb = as.numeric(as.factor(Data[PersStart, Initialization])))
  }

  ##### Call EM #####---------------------
  ## IDNames is passed to EMFunc, but EMFunc does not yet use it...not implemented yet

  out_est <- callEMFuncs(Clusters = Clusters,
                         HighestLag = HighestLag,
                         LowestLag = LowestLag,
                         Rand = Rand,
                         Rational = Rational,
                         Initialization = Initialization,
                         PreviousSol = PreviousSol,
                         IDNames = unique(Data[ , ID]),
                         K = K,
                         N = N,
                         Y = Y,
                         X = X,
                         Tni = Tni,
                         qqq = qqq,
                         nDepVar = nDepVar,
                         PersStart = PersStart,
                         PersPDiffStart = PersPDiffStart,
                         PersEnd = PersEnd,
                         PersStartU = PersStartU,
                         PersEndU = PersEndU,
                         Covariates = Covariates,
                         Conv = Conv,
                         it = it,
                         val.init = val.init,
                         smallestClN = smallestClN,
                         SigmaIncrease = SigmaIncrease,
                         call = call,
                         pbar = pbar,
                         PredictableObs = PredictableObs, 
                         Tni_NPred = Tni_NPred,
                         PersStartU_NPred = PersStartU_NPred,
                         PersEndU_NPred = PersEndU_NPred
                         )


  # ----- Parse Finish Message -----
  cat(paste0("\n LCVAR Model Estimation completed in ",  out_est$Runtime, " min"))

  # ----- Return Output Object -----

  out_est$Call$Clusters <- Clusters
  out_est$Call$Lags <- Lags

  return(out_est)


} # eof





