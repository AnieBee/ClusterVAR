
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
                        Tni,
                        qqq,
                        nDepVar,
                        PersStart,
                        PersPDiffStart,
                        PersEnd,
                        PersStartU,
                        PersEndU,
                        Covariates,
                        Conv,
                        it,
                        val.init,
                        ICType,
                        smallestClN,
                        SigmaIncrease,
                        call,
                        pbar)
{

  ### Loop over different K (# of clusters) values  -------------
  OutputAllK = vector("list", length(Clusters)) # with length(Clusters) many elements, containing the best solution for every K
  BestModels = matrix(NA, nrow = length(Clusters), ncol = 1)

  All_Solutions_for_all_starts_all_lags_all_Clusters = new.env()
  # Is in a hidden environment so it is not printed in the output unless it is specifically requested
  All_Solutions_for_all_starts_all_lags_all_Clusters$All_Solutions = vector("list", length(Clusters)) # with lenght(Clusters) many elements containing
  # the OutputListAllLags of every K, which contains all solutions for all Lags and all starts for that K


  ### Names for Output ### ######
  OutputAllKnames = NULL
  for (K in Clusters)
  {
    OutputAllKnames = c(OutputAllKnames,
                        paste(c("Best solution of all models where the number of clusters is:", K), collapse = " "))
  }

  names(OutputAllK) = OutputAllKnames
  #######


  # ----- Compute Lag Combination for each Cluster (Jonas) -----
  # We do this already up here, so we can use it to normalize the progress bar
  l_LagsPerCluster <- list()
  for(K in Clusters)  l_LagsPerCluster[[K]] <- calculateLagList(K = K,
                                                                HighestLag = HighestLag,
                                                                LowestLag = LowestLag)
  # Compute number of fitted models
  n_Fit <- sum(unlist(lapply(l_LagsPerCluster, nrow))) *  (Rand + as.numeric(Rational) + as.numeric(!is.null(Initialization)) + as.numeric(PreviousSol)) # number of models times number of random restarts

  # browser()

  # ----- Create Progress Bar (Jonas) -----
  # Calculate Max

  # Set up progress bar
  if(pbar==TRUE) pb <- txtProgressBar(min = 0, max=n_Fit, initial=0, char="-", style = 3)
  pb_counter <- 0
  start_time <- proc.time()[3]

  ClustCount = 1
  for(K in Clusters)  {

    # cat(c("\n", K, "Clusters: "))
    # cat(c("\n", K, "Clusters: "), file = "EMwarnings.txt", append = TRUE)

    # LagsList = calculateLagList(K = K,
    #                             HighestLag =
    #                               HighestLag,
    #                             LowestLag = LowestLag)

    # if(K == 1) LagsList = t(LagsList)  # dimension changes when number of clusters is equal to 1 >> In that case, transpose to get same dimension again

    LagCombinations <- nrow(l_LagsPerCluster[[K]]) #dim(LagsList)[1]
    LagsList <- l_LagsPerCluster[[K]] # Jonas: now already computed above

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


    # Jonas: ???
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
                                                                                 PersPDiffStart = PersPDiffStart)


    # browser()

    EMCallVec = c(as.character(1:Rand),
                  ifelse(Rational, "Rational", NULL),
                  switch(!is.null(Initialization), "FALSE" = "Initialization"), # default is NULL
                  ifelse(PreviousSol, "Previous", NULL))

    usePrevLagSol = FALSE   # Make sure the previous lag solution-using start is not called before a previous solution exists
    for(LagCounter in 1:LagCombinations) {

      # cat(c('\n', " Lags:", LagsList[LagCounter, ], '\n'))
      # cat(c('\n', " Lags:", LagsList[LagCounter, ], '\n'), file = "EMwarnings.txt", append = TRUE)
      ### Initialization Prerequesites: calcuate coefficeints passed to initial clustering solutions ###---------------------------------
      #PersPDiffStart, PersStartU etc can (must) be integers instead of vectors in calculateCoefficientsForRandoAndRational


      #### Random starts ###
      StartCounter = 0
      while (StartCounter != length(EMCallVec))
        # ToDo: in all EMFunc IDNames is passed but not used
      {

        StartCounter = StartCounter + 1

        # cat(c("*", EMCallVec[StartCounter], "  "))
        # cat(c("*", EMCallVec[StartCounter], "  "), file = "EMwarnings.txt", append = TRUE)

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
                                                       [[PrevBestRun[1]]][[PrevBestRun[2]]]$Classification,
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

                               Y = Y, X = X, Lags = LagsList[LagCounter, ], K = K, N = N, Tni = Tni, qqq = qqq, nDepVar = nDepVar,
                               PersStart = PersStart, PersPDiffStart = PersPDiffStart, PersEnd = PersEnd,
                               PersStartU = PersStartU, PersEndU = PersEndU,
                               Covariates = Covariates, smallestClN = smallestClN, SigmaIncrease = SigmaIncrease),
                 IDNames = IDNames, Y = Y, X = X, K = K, N = N, Tni = Tni, qqq = qqq, nDepVar = nDepVar,
                 PersPDiffStart = PersPDiffStart, PersEnd = PersEnd, PersStartU = PersStartU, PersEndU = PersEndU,
                 Covariates = Covariates, Conv = Conv, it = it, smallestClN = smallestClN, ICType = ICType,
                 SigmaIncrease = SigmaIncrease)

        FitAllLags[LagCounter, StartCounter] = OutputListAllLags[[LagCounter]][[StartCounter]]$IC
        
        # ----- Update progress bar (Jonas) -----
        pb_counter <- pb_counter + 1 # this should be the inner loop, so that should be what I need
        if(pbar==TRUE) setTxtProgressBar(pb, pb_counter)

        ########## Addition to speed up estimation in case K == 1
        #If K = 1 all partitions of individuals are the same (all people in the same cluster)
        # so all starts will lead to the exact same outcome
        if (K == 1) {
          while (StartCounter != length(EMCallVec))
            # ToDo: in all EMFunc IDNames is passed but not used
          {
            StartCounter = StartCounter + 1

            # cat(c("*", EMCallVec[StartCounter], "  "))
            # cat(c("*", EMCallVec[StartCounter], "  "), file = "EMwarnings.txt", append = TRUE)

            OutputListAllLags[[LagCounter]][[StartCounter]] = OutputListAllLags[[LagCounter]][[StartCounter - 1]]
            FitAllLags[LagCounter, StartCounter] = FitAllLags[LagCounter, StartCounter - 1]
            
            # ----- Update progress bar (Jonas) -----
            pb_counter <- pb_counter + 1 # this should be the inner loop, so that should be what I need
            if(pbar==TRUE) setTxtProgressBar(pb, pb_counter)
          }
        }
        ################################################################

      } # End of Start loop

      PrevBestRun = arrayInd(which.min(FitAllLags), dim(FitAllLags))
      usePrevLagSol = TRUE
      LagCounter = LagCounter + 1

    } # End of Lag loop

    BestRunOneK = arrayInd(which.min(FitAllLags), dim(FitAllLags))

    # ToDo: update Classification with ID
    All_Solutions_for_all_starts_all_lags_all_Clusters$All_Solutions[[ClustCount]] = OutputListAllLags
    ModelCall = list(Clusters = Clusters, Lags = LowestLag:HighestLag,
                     Rand = Rand, Rational = Rational, Initialization = Initialization,
                     ICType = ICType, Covariates = Covariates)
    OutputAllK[[ClustCount]] = OutputListAllLags[[BestRunOneK[1]]][[BestRunOneK[2]]]
    BestModels[ClustCount, ] =  paste(c("According to the selected infromation criterion, the best solution of all models where the number of clusters is", K,
                                        "has a lag order of:", OutputListAllLags[[BestRunOneK[1]]][[BestRunOneK[2]]]$Lags), collapse = " ")
    ClustCount = ClustCount + 1

  } # End of K cluster loop


  # ----- Timer -----
  end_time <- proc.time()[3] - start_time
  end_time_min <- round(end_time/60, 2)


  outlist <- list(Best_Models = BestModels,
                  Best_solutions_per_cluster = OutputAllK,
                  call = call,
                  Call = ModelCall,
                  All_Solutions_for_all_starts_all_lags_all_clusters = All_Solutions_for_all_starts_all_lags_all_Clusters,
                  Runtime = end_time_min)

  class(outlist) <- c("ClusterVAR", class(outlist))

  return(outlist)

}  # eoF








