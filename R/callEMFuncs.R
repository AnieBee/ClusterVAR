
callEMFuncs <- function(Clusters,
                        HighestLag,
                        LowestLag,
                        Rand,
                        Rational,
                        Initialization,
                        PreviousSol,
                        IDNames,
                        K,
                        N,
                        Y,
                        X,
                        qqq,
                        nDepVar,
                        PersStart,
                        PersEnd,
                        Covariates,
                        Conv,
                        it,
                        val.init,
                        smallestClN,
                        SigmaIncrease,
                        pbar,
                        NewPredictableObs,
                        LaggedPredictObs,
                        PredictableObsConc,
                        LaggedPredictObsConc,
                        Tni_NPred,
                        PersStartU_NPred,
                        PersEndU_NPred)
{

  ### Loop over different K (# of clusters) values  -------------
  All_Solutions = vector("list", length(Clusters)) # with lenght(Clusters) many elements containing
  # the OutputListAllLags of every K, which contains all solutions for all Lags and all starts for that K


  # ----- Compute Lag Combination for each Cluster -----
  # We do this already up here, so we can use it to normalize the progress bar
  l_LagsPerCluster <- list()
  for(K in Clusters)  l_LagsPerCluster[[K]] <- calculateLagList(K = K,
                                                                HighestLag = HighestLag,
                                                                LowestLag = LowestLag)
  # Compute number of fitted models
  n_Fit <- sum(unlist(lapply(l_LagsPerCluster, nrow))) * (Rand + as.numeric(Rational) + as.numeric(!is.null(Initialization)) + as.numeric(PreviousSol)) # number of models times number of random restarts



  # ----- Create Progress Bar -----
  # Set up progress bar
  if(pbar==TRUE) pb <- txtProgressBar(min = 0, max=n_Fit, initial=0, char="-", style = 3)
  pb_counter <- 0
  start_time <- proc.time()[3]

  ClustCount = 1
  for(K in Clusters)  {

    LagCombinations <- nrow(l_LagsPerCluster[[K]]) #dim(LagsList)[1]
    LagsList <- l_LagsPerCluster[[K]]

    # OutputListAllLags[[Lags]][[Start]]
    OutputListAllLags = createOutputList(LagCombinations = LagCombinations,
                                         Rand = Rand,
                                         Rational = Rational,
                                         Initialization = Initialization,
                                         PreviousSol = PreviousSol)
    # Fit[Lags, Start]
    FitAllLags = array(NA, dim = c(LagCombinations, Rand + as.numeric(Rational) +
                                     as.numeric(!is.null(Initialization))
                                   + as.numeric(PreviousSol))) # to store fit of output


    ### Initialization Prerequesites: calcuate coefficeints passed to initial clustering solutions ###---------------------------------
    CoeffsForRandoAndRationalList = callCalculateCoefficientsForRandoAndRational(Covariates = Covariates,
                                                                                 K = K,
                                                                                 N = N,
                                                                                 nDepVar = nDepVar,
                                                                                 qqq = qqq,
                                                                                 HighestLag = HighestLag,
                                                                                 LowestLag = LowestLag,
                                                                                 PersEnd = PersEnd,
                                                                                 PersStart = PersStart,
                                                                                 Y = Y,
                                                                                 X = X,
                                                                                 NewPredictableObs = NewPredictableObs,
                                                                                 LaggedPredictObs = LaggedPredictObs)
    #EMCallVec[StartCounter] gives name/type of current start (e.g., random start nr "4", or "Previous")
    EMCallVec = c(as.character(1:Rand),
                  ifelse(Rational, "Rational", NULL),
                  switch(!is.null(Initialization), "FALSE" = "Initialization"), # default is NULL
                  ifelse(PreviousSol, "Previous", NULL))

    usePrevLagSol = FALSE   # Make sure the previous lag solution-using start is not called before a previous solution exists
    for(LagCounter in 1:LagCombinations) {


         #### Random starts ###
      StartCounter = 0
      while (StartCounter != length(EMCallVec))
        # ToDo: in all EMFunc IDNames is passed but not used
      {

        StartCounter = StartCounter + 1

        OutputListAllLags[[LagCounter]][[StartCounter]] =
            EMFunc(Init = EMInit(InitMT =
                                 switch(EMCallVec[StartCounter],
                                        "Rational" = InitRat(K = K,
                                                             CoefficientsForRandoAndRational =
                                                               CoeffsForRandoAndRationalList[[max(LagsList[LagCounter, ])]]),

                                        "Initialization" = val.init,

                                        "Previous" = if(usePrevLagSol)
                                        {   ## Previous Sol ##
                                          t(dummy_cols(OutputListAllLags
                                                       [[PrevBestRun[1]]][[PrevBestRun[2]]]$Classification[1, ],
                                                       remove_first_dummy = FALSE)[ , -1])
                                        }
                                        else
                                        {   ## PseudoRand (a previous sol does not exist yet) ##
                                          InitPseudoRand(N = N, K = K, smallestClN = smallestClN,
                                                         CoefficientsForRandoAndRational =
                                                           CoeffsForRandoAndRationalList[[max(LagsList[LagCounter, ])]])
                                        },
                                        ### Default: PseudoRand initialization ###
                                        InitPseudoRand(N = N, K = K, smallestClN = smallestClN,
                                                       CoefficientsForRandoAndRational =
                                                         CoeffsForRandoAndRationalList[[max(LagsList[LagCounter, ])]])
                                 ),

                               Y = Y, X = X, Lags = LagsList[LagCounter, ], K = K, N = N, qqq = qqq, nDepVar = nDepVar,
                               PersStart = PersStart, PersEnd = PersEnd,
                               NewPredictableObs = NewPredictableObs,
                               LaggedPredictObs = LaggedPredictObs,
                               PredictableObsConc = PredictableObsConc,
                               LaggedPredictObsConc = LaggedPredictObsConc,
                               Tni_NPred = Tni_NPred,
                               PersStartU_NPred = PersStartU_NPred,
                               PersEndU_NPred = PersEndU_NPred,
                               Covariates = Covariates, smallestClN = smallestClN, SigmaIncrease = SigmaIncrease),
                 IDNames = IDNames, Y = Y, X = X, K = K, N = N, Tni_NPred = Tni_NPred, qqq = qqq, nDepVar = nDepVar,
                 NewPredictableObs = NewPredictableObs, LaggedPredictObs = LaggedPredictObs,
                 PredictableObsConc = PredictableObsConc,
                 LaggedPredictObsConc = LaggedPredictObsConc,
                 PersEnd = PersEnd, PersStart = PersStart,
                 PersStartU_NPred = PersStartU_NPred, PersEndU_NPred = PersEndU_NPred,
                 Covariates = Covariates, Conv = Conv, it = it, smallestClN = smallestClN,
                 SigmaIncrease = SigmaIncrease)

        FitAllLags[LagCounter, StartCounter] = OutputListAllLags[[LagCounter]][[StartCounter]]$SC

        # ----- Update progress bar  -----
        pb_counter <- pb_counter + 1
        if(pbar==TRUE) setTxtProgressBar(pb, pb_counter)

        ########## Addition to speed up estimation in case K == 1
        #If K = 1 all partitions of individuals are the same (all people in the same cluster)
        # so all starts will lead to the exact same outcome
        if (K == 1) {
          while (StartCounter != length(EMCallVec))
            # ToDo: in all EMFunc IDNames is passed but not used
          {
            StartCounter = StartCounter + 1

            OutputListAllLags[[LagCounter]][[StartCounter]] = OutputListAllLags[[LagCounter]][[StartCounter - 1]]
            FitAllLags[LagCounter, StartCounter] = FitAllLags[LagCounter, StartCounter - 1]

            # ----- Update progress bar  -----
            pb_counter <- pb_counter + 1
            if(pbar==TRUE) setTxtProgressBar(pb, pb_counter)
          }
        }
      } # End of Start loop

      PrevBestRun = arrayInd(which.min(FitAllLags), dim(FitAllLags))
      usePrevLagSol = TRUE

    } # End of Lag loop

    All_Solutions[[ClustCount]] = OutputListAllLags
    ModelCall = list(Clusters = Clusters, Lags = LowestLag:HighestLag,
                     Rand = Rand, Rational = Rational, Initialization = Initialization,
                     Covariates = Covariates)
    ClustCount = ClustCount + 1

  } # End of K cluster loop


  # ----- Timer -----
  end_time <- proc.time()[3] - start_time
  end_time_min <- round(end_time/60, 2)


  outlist <- list(Call = ModelCall,
                  All_Models = All_Solutions,
                  Runtime = end_time_min)

  class(outlist) <- c("ClusterVAR", class(outlist))

  return(outlist)

}  # eoF








