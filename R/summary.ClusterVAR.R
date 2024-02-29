
summary.ClusterVAR <- function(object,
                               ...) {

  # Fill in defaults
  args <- list(...)
  if(is.null(args$show)) show <- "BPC" else show <- args$show
  if(is.null(args$TS_criterion)) TS_criterion <- "SC" else TS_criterion <- args$TS_criterion
  if(is.null(args$global_criterion)) global_criterion <- "BIC" else global_criterion <- args$global_criterion
  if(is.null(args$Number_of_Clusters)) Number_of_Clusters <- NULL else Number_of_Clusters <- args$Number_of_Clusters

  #show =  c("BPC", "BO", "GNC") = c("Best-per-number-of-clusters", "Best-overall", "Given-a-number-of-clusters")
  #TS_criterion = c("SC", "HQ")
  #Global_criterion = c("BIC", "ICL")
  # Number_of_Clusters is only used if show == "GNC"

  # -------------Check all input is as expected---------------
  if (!(show %in% c("BPC", "BO", "GNC"))) {
    stop("Invalid value for 'show'. Please choose from: ", paste(c("BPC", "BO", "GNC"), collapse = ", "))
  }
  if (!(TS_criterion %in% c("SC", "HQ"))) {
    stop("Invalid value for 'TS_criterion'. Please choose from: ", paste(c("SC", "HQ"), collapse = ", "))
  }
  if (!(global_criterion %in% c("BIC", "ICL"))) {
    stop("Invalid value for 'global_criterion'. Please choose from: ", paste(c("BIC", "ICL"), collapse = ", "))
  }
  if (show == "GNC" && is.null(Number_of_Clusters)) {
    stop("If 'show' is 'GNC', you must specify a value for 'Number_of_Clusters'.")
  }
  if(show == "GNC" && !(Number_of_Clusters %in% object$Call$Clusters)){
    stop("The value you specified for 'Number_of_Clusters' is not contained in the Cluster Sequence of your fitted LCVAR Model Object.")
  }
  # -------------

  l_LagsPerCluster <- list()
  for(K in object$Call$Clusters)  l_LagsPerCluster[[K]] <- calculateLagList(K = K, HighestLag = max(object$Call$Lags), LowestLag = min(object$Call$Lags))
  NumberStarts = (object$Call$Rand + as.numeric(object$Call$Rational) + as.numeric(!is.null(object$Call$Initialization)) + 1)
  # last element is for PreviousSol which is TRUE as a default

  # -------------Create one of the three different types of FunctionOutput---------------

  if(show == "GNC"){

    clusterCounter = which(object$Call$Clusters == Number_of_Clusters)
    GivenClusterOutput = object$Models[[clusterCounter]]
    LagCombinations = nrow(l_LagsPerCluster[[Number_of_Clusters]])

    # Fit[Lags, Start]
    FitAllLags = array(NA, dim = c(LagCombinations, NumberStarts)) # to store fit of output
    FunctionOutput = data.frame(matrix(NA, nrow = LagCombinations, ncol = 8),
                                row.names = apply(l_LagsPerCluster[[Number_of_Clusters]], 1, function(x) paste(c("Lags", x), collapse = " ")))
    colnames(FunctionOutput) = c("log-likelihood", "parameters", "SC", "HQ", "BIC", "ICL", "Converged", "Proportions")


    for(LagCounter in 1:LagCombinations) {
      for(StartCounter in 1:NumberStarts){
        FitAllLags[LagCounter, StartCounter] = switch(TS_criterion,
                                                      "SC" = GivenClusterOutput[[LagCounter]][[StartCounter]]$SC,
                                                      "HQ" = GivenClusterOutput[[LagCounter]][[StartCounter]]$HQ)

      }
      BestModel = which.min(FitAllLags[LagCounter, ]) # Best model for this LagCounter across all starts
      #ModelNamesLags[LagCounter] = paste(GivenClusterOutput[[LagCounter]][[BestModel]]$Lags, collapse = " ")
      FunctionOutput[LagCounter, "log-likelihood"] = GivenClusterOutput[[LagCounter]][[BestModel]]$last.loglik
      FunctionOutput[LagCounter, "parameters"] = GivenClusterOutput[[LagCounter]][[BestModel]]$nParameters
      FunctionOutput[LagCounter, "SC"] = GivenClusterOutput[[LagCounter]][[BestModel]]$SC
      FunctionOutput[LagCounter, "HQ"] = GivenClusterOutput[[LagCounter]][[BestModel]]$HQ
      FunctionOutput[LagCounter, "BIC"] = GivenClusterOutput[[LagCounter]][[BestModel]]$BIC
      FunctionOutput[LagCounter, "ICL"] = GivenClusterOutput[[LagCounter]][[BestModel]]$ICL
      FunctionOutput[LagCounter, "Converged"] = GivenClusterOutput[[LagCounter]][[BestModel]]$Converged
      FunctionOutput[LagCounter, "Proportions"] = paste(round(GivenClusterOutput[[LagCounter]][[BestModel]]$Proportions, 2), collapse = " ")
    }
    message = paste0(c("---------------------------------------------------\n",
                       "All lags for number of clusters =", Number_of_Clusters,
                       "\n---------------------------------------------------\n"))
  } else {
    # The below calculation calculates the solution for (show == "BPC") but all calculation steps are also needed if (show == "BO")
    # For each number of of clusters, find the single best-fitting time-series model for each cluster number

    FunctionOutput = data.frame(matrix(NA, nrow = length(object$Call$Clusters), ncol = 7),
                                row.names = apply(as.matrix(object$Call$Clusters), 1, function(x) paste(c(x, "Clusters"), collapse = " ")))
    colnames(FunctionOutput) = c(paste(c("Lags selected by", TS_criterion), collapse = " "),
                                 "log-likelihood", "parameters", "BIC", "ICL", "Converged", "Proportions")

    # browser()

    ClustCount = 1
    for(K in object$Call$Clusters){

      LagCombinations = nrow(l_LagsPerCluster[[K]])
      FitAllLags = array(NA, dim = c(LagCombinations, NumberStarts)) # to store fit of output

      for(LagCounter in 1:LagCombinations) {
        for(StartCounter in 1:NumberStarts){
          FitAllLags[LagCounter, StartCounter] = switch(TS_criterion,
                                                        "SC" = object$Models[[ClustCount]][[LagCounter]][[StartCounter]]$SC,
                                                        "HQ" = object$Models[[ClustCount]][[LagCounter]][[StartCounter]]$HQ)
        }
      }
      BestRunOneK = arrayInd(which.min(FitAllLags), dim(FitAllLags))
      FunctionOutput[ClustCount, 1] = paste(object$Models[[ClustCount]][[BestRunOneK[1]]][[BestRunOneK[2]]]$Lags, collapse = " ")
      FunctionOutput[ClustCount, "log-likelihood"] = object$Models[[ClustCount]][[BestRunOneK[1]]][[BestRunOneK[2]]]$last.loglik
      FunctionOutput[ClustCount, "parameters"] = object$Models[[ClustCount]][[BestRunOneK[1]]][[BestRunOneK[2]]]$nParameters
      FunctionOutput[ClustCount, "BIC"] = object$Models[[ClustCount]][[BestRunOneK[1]]][[BestRunOneK[2]]]$BIC
      FunctionOutput[ClustCount, "ICL"] = object$Models[[ClustCount]][[BestRunOneK[1]]][[BestRunOneK[2]]]$ICL
      FunctionOutput[ClustCount, "Converged"] = object$Models[[ClustCount]][[BestRunOneK[1]]][[BestRunOneK[2]]]$Converged
      FunctionOutput[ClustCount, "Proportions"] = paste(round(object$Models[[ClustCount]][[BestRunOneK[1]]][[BestRunOneK[2]]]$Proportions, 2), collapse = " ")


      ClustCount = ClustCount + 1
    }


    if(show == "BPC"){
      message = paste0(c("---------------------------------------------------\n",
                         "The best lags for any number of clusters as selected by the", TS_criterion,
                         "\n---------------------------------------------------\n"))
      # FunctionOutput stays as calculated above
    }
    if (show == "BO"){
      BestOverall = switch(global_criterion,
                           "BIC" = which.min(FunctionOutput$BIC),
                           "ICL" = which.min(FunctionOutput$ICL))
      message = paste0(c("---------------------------------------------------\n",
                         "The best lags for any number of clusters as selected by the", TS_criterion,
                         "\n",
                         "The best number of clusters as selected by the", global_criterion,
                         "\n---------------------------------------------------\n"))
      FunctionOutput = FunctionOutput[BestOverall, ]
    }

  } # end if: what to show?
  #### Return FunctionOutput here ##
  FunctionOutput = list(message = message, FunctionOutput = FunctionOutput)
  class(FunctionOutput) <- c("ClusterVARSummary", class(FunctionOutput))

  return(FunctionOutput)


} # eoF
