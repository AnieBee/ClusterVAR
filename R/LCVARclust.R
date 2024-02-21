

LCVARclust <- function(Data,
                       yVars,
                       Beep,
                       Day = NULL,
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
# If each measurement is done on the same day, don't specify day but only specify beep.
# Only if the first measurment on each day should be removed from calculations, the 
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
  stopifnot(HighestLag < 3) # highest lag number that is allowed is 3
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

    PreviousSol = TRUE # Use solution of previous Lags as a start

  # Y has to be numeric, X can be numeric or factor
  # Y is data of form m \times (sum^N (nObs) )
  # ID = has to be factor
  ##### Preprocessing of Data Set #####--------------------
  Data = as.data.frame(Data)
    if(is.null(Day)){
        Data = Data[order(Data[ , ID], Data[ , Beep]), ]
    }else{
        Data = Data[order(Data[ , ID], Data[ , Day], Data[ , Beep]), ]
    }
  # order Data according to ID, make sure an individual's observations occur one after another with first obs first, second second etc
  # observations have to occur ascending in time

  # Remove rows with NA values
  Data <- Data[complete.cases(Data), ]
  
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
  stopifnot(is.numeric(Y))
  stopifnot(N >= (smallestClN * max(Clusters)))# check that for the highest number of clusters requested
  # there can be at least the smallest permissable number of people per cluster

  # Whole Sample of Obs
  PersStart = cumsum(c(0, nObs[-length(nObs)])) + 1 # Start of individual time series for every pers in Y or W
  PersEnd = cumsum(nObs) # end of individual time series of every pers, last obs of every pers in Y or W

  
  ### ------------ whether an observation in Y can be predicted (PredictableObs) ------------------
  PredictableObs = vector("list", HighestLag)  # from 1:HighestLag but only LowestLag:HighestLag elements are filled
  for (lagRunner in LowestLag:HighestLag) 
  {
      ind_lag_all <- lapply(1:N, function(x) {
          dayArgument_pers = if (!is.null(Day)) Data[PersStart[x]:PersEnd[x], Day] else NULL
          out <- DetPredSubj(beep = Data[PersStart[x]:PersEnd[x], Beep], 
                             day = dayArgument_pers, 
                             MaxLag = lagRunner)
          return(out)
      })
      ind_lag_all <- do.call(c, ind_lag_all)
      stopifnot(length(which(is.na(ind_lag_all))) == 0) # check that everything went well
      PredictableObs[[lagRunner]] = rbind(ind_lag_all, Pers)
  }
 ###### End of determining PredictableObs ------------------------------
  
  
  ### Runners that depend on the number of predictable obs (and thus on the number of lags) ### ----
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
          # Tni_NPred[[lagRunner]] # length of U per person,
          # With presample of first Lags-obs removed for every individual and all other obs that cannot be predicted removed,
        }
      ## Check each person has at least 5 predictable observations
      stopifnot(all(Tni_NPred[[lagRunner]] > 4))
      
      ## Change PredictableObs so that it contains the position of the observations in Y, and W, and X, that can be predicted
      PredictableObs[[lagRunner]] = which(PredictableObs[[lagRunner]][1,] ==  1)
     
      PersStartU_NPred[[lagRunner]] = cumsum(c(0, Tni_NPred[[lagRunner]][-length(Tni_NPred[[lagRunner]])])) + 1 
      # Start of individual time series for every pers in U ( = non-predictables are removed)
      PersEndU_NPred[[lagRunner]] = cumsum(Tni_NPred[[lagRunner]]) # Last obs of a pers in U
      
  }
  
  ### Create new PredictableObs
  NewPredictableObs = vector("list", HighestLag)  # from 1:HighestLag but only LowestLag:HighestLag elements are filled
  LaggedPredictObs = vector("list", HighestLag)  # from 1:HighestLag but only LowestLag:HighestLag elements are filled
  GetSequences <- function(trunner, Lag) {
      lapply(trunner, function(t) (t - 1):(t - Lag))
  }
  for (lagRunner in LowestLag:HighestLag)
  {
      NewPredictableObsSmall = vector("list", N)
      LaggedPredictObsSmall = vector("list", N)
      for (i in 1:N) {
          NewPredictableObsSmall[[i]] = c(intersect(PredictableObs[[lagRunner]], c(PersStart[i]:PersEnd[i])))
          LaggedPredictObsSmall[[i]]  = unlist(GetSequences(NewPredictableObsSmall[[i]], lagRunner))
      }
      NewPredictableObs[[lagRunner]] = NewPredictableObsSmall
      LaggedPredictObs[[lagRunner]] = LaggedPredictObsSmall
  }
  
##########################
 
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
                         qqq = qqq,
                         nDepVar = nDepVar,
                         PersStart = PersStart,
                         PersEnd = PersEnd,
                         Covariates = Covariates,
                         Conv = Conv,
                         it = it,
                         val.init = val.init,
                         smallestClN = smallestClN,
                         SigmaIncrease = SigmaIncrease,
                         call = call,
                         pbar = pbar,
                         PredictableObs = PredictableObs,
                         NewPredictableObs = NewPredictableObs, 
                         LaggedPredictObs = LaggedPredictObs,
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





