
summary.ClusterVAR <- function(object, show = "BPC",  TS_criterion = "SC",
                               global_criterion = "BIC",
                               Number_of_Clusters = NULL,
                               Number_of_Lags = NULL,
                               ...) {

  # Fill in defaults
  # args <- list(...)
  # if(is.null(args$show)) show <- "BPC" else show <- args$show
  # if(is.null(args$TS_criterion)) TS_criterion <- "SC" else TS_criterion <- args$TS_criterion
  # if(is.null(args$global_criterion)) global_criterion <- "BIC" else global_criterion <- args$global_criterion
  # if(is.null(args$Number_of_Clusters)) Number_of_Clusters <- NULL else Number_of_Clusters <- args$Number_of_Clusters
  # if(is.null(args$Number_of_Lags)) Number_of_Lags <- min(object$Call$Lags) else Number_of_Lags <- args$Number_of_Lags
  if(is.null(Number_of_Lags)) Number_of_Lags <- min(object$Call$Lags)

  #TS_criterion = c("SC", "HQ")
  #Global_criterion = c("BIC", "ICL")
  # Number_of_Clusters is only used if show == "GNC"

  # -------------Check all input is as expected---------------
  if (!(show %in% c("BPC", "GNL", "GNC"))) {
    stop("Invalid value for 'show'. Please choose from: ", paste(c("BPC", "GNL", "GNC"), collapse = ", "))
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
    stop("The value you specified for 'Number_of_Clusters' is not contained in the Cluster sequence of your fitted LCVAR Model Object.")
  }
  if(show == "GNL" && !(Number_of_Lags %in% object$Call$Lags)){
    stop("The value you specified for 'Number_of_Lags' is not contained in the Lag sequence of your fitted LCVAR Model Object.")
  }
  # -------------

  l_LagsPerCluster <- list()
  for(K in object$Call$Clusters)  l_LagsPerCluster[[K]] <- calculateLagList(K = K, HighestLag = max(object$Call$Lags), LowestLag = min(object$Call$Lags))
  NumberStarts = (object$Call$Rand + as.numeric(object$Call$Rational) + as.numeric(!is.null(object$Call$Initialization)) + 1)
  # last element is for PreviousSol which is TRUE as a default

  # -------------Create one of the three different types of FunctionOutput---------------

  if(show == "GNC"){

    LagCombinations = nrow(l_LagsPerCluster[[Number_of_Clusters]])
    FunctionOutput = data.frame(matrix(NA, nrow = LagCombinations, ncol = 7),
                                row.names = apply(l_LagsPerCluster[[Number_of_Clusters]], 1, function(x) paste(c("Lags", x), collapse = " ")))
    colnames(FunctionOutput) = c("log-likelihood", "parameters", "Lags", "SC", "HQ", "Converged", "Proportions")


    for(LagCounter in 1:LagCombinations) {
      BestModel = coef(object, Model = l_LagsPerCluster[[Number_of_Clusters]][LagCounter, ])  # Best model for this LagCounter across all starts (based on likelihood)
      ExtractedLags = BestModel$Lags
      OrderedLags = ExtractedLags[order(ExtractedLags)]
      OrderedProportions = BestModel$Proportions[order(ExtractedLags)]
      FunctionOutput[LagCounter, "log-likelihood"] = BestModel$last.loglik
      FunctionOutput[LagCounter, "parameters"] = BestModel$nParameters
      FunctionOutput[LagCounter, "Lags"] = paste(OrderedLags, collapse = " ")
      FunctionOutput[LagCounter, "SC"] = BestModel$SC
      FunctionOutput[LagCounter, "HQ"] = BestModel$HQ
      FunctionOutput[LagCounter, "Converged"] = BestModel$Converged
      FunctionOutput[LagCounter, "Proportions"] = paste(round(OrderedProportions, 2), collapse = " ")
    }
    BestOverall = switch(TS_criterion,
                         "SC" = which.min(FunctionOutput$SC),
                         "HQ" = which.min(FunctionOutput$HQ))
    message = paste0(c("---------------------------------------------------\n",
                       "All lags for number of clusters =", Number_of_Clusters,
                       "\n",
                       "For this number of clusters the", TS_criterion,"selects:", row.names(FunctionOutput)[BestOverall],
                       "\n---------------------------------------------------\n"))

  }
  if(show == "BPC"){
    # The below calculates the solution for (show == "BPC")
    # For each number of clusters, find the single best-fitting time-series model for each cluster number

    FunctionOutput = data.frame(matrix(NA, nrow = length(object$Call$Clusters), ncol = 7),
                                row.names = apply(as.matrix(object$Call$Clusters), 1, function(x) paste(c(x, "Clusters"), collapse = " ")))
    colnames(FunctionOutput) = c(paste(c("Lags selected by", TS_criterion), collapse = " "),
                                 "log-likelihood", "parameters", "BIC", "ICL", "Converged", "Proportions")

    ClustCount = 1
    for(K in object$Call$Clusters){

      LagCombinations = nrow(l_LagsPerCluster[[K]])
      FitAllLags = array(NA, dim = c(LagCombinations, 2)) # to store fit of output
      FitStartsWithinLag = array(NA, dim = c(NumberStarts))

      for(LagCounter in 1:LagCombinations) {
        for(StartCounter in 1:NumberStarts){
          FitStartsWithinLag[StartCounter] = object$All_Models[[ClustCount]][[LagCounter]][[StartCounter]]$last.loglik
        }
        FitAllLags[LagCounter, 1] = which.max(FitStartsWithinLag)[1] # determines the best start for each lag based on likelihood
        FitAllLags[LagCounter, 2] = switch(TS_criterion,
                                                  "SC" = object$All_Models[[ClustCount]][[LagCounter]][[FitAllLags[LagCounter, 1]]]$SC,
                                                  "HQ" = object$All_Models[[ClustCount]][[LagCounter]][[FitAllLags[LagCounter, 1]]]$HQ)
      }
      BestRunOneK = which.min(FitAllLags[, 2])[1]
      BPCFinalModel = object$All_Models[[ClustCount]][[BestRunOneK]][[FitAllLags[BestRunOneK, 1]]]
      ExtractedLags = BPCFinalModel$Lags
      OrderedLags = ExtractedLags[order(ExtractedLags)]
      OrderedProportions = BPCFinalModel$Proportions[order(ExtractedLags)]
      FunctionOutput[ClustCount, 1] = paste(OrderedLags, collapse = " ")
      FunctionOutput[ClustCount, "log-likelihood"] = BPCFinalModel$last.loglik
      FunctionOutput[ClustCount, "parameters"] = BPCFinalModel$nParameters
      FunctionOutput[ClustCount, "BIC"] = BPCFinalModel$BIC
      FunctionOutput[ClustCount, "ICL"] = BPCFinalModel$ICL
      FunctionOutput[ClustCount, "Converged"] = BPCFinalModel$Converged
      FunctionOutput[ClustCount, "Proportions"] = paste(round(OrderedProportions, 2), collapse = " ")

      ClustCount = ClustCount + 1

    }

    # remove BIC and ICL from this output, because models with unequal number of lags cannot be compared to one another
    FunctionOutput = FunctionOutput[, -c(2, 4, 5)] # remove liklike, BIC and ICL because they are not comparable for all models
    message = paste0(c("---------------------------------------------------\n",
                       "The best lags for each number of clusters as selected by the", TS_criterion,
                       "\n---------------------------------------------------\n"))
    # FunctionOutput stays as calculated above
  }
  if (show == "GNL"){

    FunctionOutput = data.frame(matrix(NA, nrow = length(object$Call$Clusters), ncol = 7),
                                row.names = apply(as.matrix(object$Call$Clusters), 1, function(x) paste(c(x, "Clusters"), collapse = " ")))
    colnames(FunctionOutput) = c(paste(c("Lags"), collapse = " "),
                                 "log-likelihood", "parameters", "BIC", "ICL", "Converged", "Proportions")

    ClustCount = 1
    for(K in object$Call$Clusters){

      LagCounterresult <- apply(l_LagsPerCluster[[K]], 1, function(row) all(row == Number_of_Lags))
      LagCounter <- which(LagCounterresult)

      FitAllLags = array(NA, dim = c(NumberStarts)) # to store fit of output

      for(StartCounter in 1:NumberStarts){
        FitAllLags[StartCounter] = object$All_Models[[ClustCount]][[LagCounter]][[StartCounter]]$last.loglik
      }
      BestRunOneK = which.max(FitAllLags)[1]
      ExtractedLags = object$All_Models[[ClustCount]][[LagCounter]][[BestRunOneK]]$Lags # Don't need to be ordered because they are all ordered by default of being all the same lag
      FunctionOutput[ClustCount, 1] = paste(ExtractedLags, collapse = " ")
      FunctionOutput[ClustCount, "log-likelihood"] = object$All_Models[[ClustCount]][[LagCounter]][[BestRunOneK]]$last.loglik
      FunctionOutput[ClustCount, "parameters"] = object$All_Models[[ClustCount]][[LagCounter]][[BestRunOneK]]$nParameters
      FunctionOutput[ClustCount, "BIC"] = object$All_Models[[ClustCount]][[LagCounter]][[BestRunOneK]]$BIC
      FunctionOutput[ClustCount, "ICL"] = object$All_Models[[ClustCount]][[LagCounter]][[BestRunOneK]]$ICL
      FunctionOutput[ClustCount, "Converged"] = object$All_Models[[ClustCount]][[LagCounter]][[BestRunOneK]]$Converged
      FunctionOutput[ClustCount, "Proportions"] = paste(round(object$All_Models[[ClustCount]][[LagCounter]][[BestRunOneK]]$Proportions, 2), collapse = " ")

      ClustCount = ClustCount + 1
    }

    BestOverall = switch(global_criterion,
                         "BIC" = which.min(FunctionOutput$BIC),
                         "ICL" = which.min(FunctionOutput$ICL))
    message = paste0(c("---------------------------------------------------\n",
                       "All models where all lags =", Number_of_Lags,
                       "\n",
                       "For this number of lags the", global_criterion,"selects:", row.names(FunctionOutput)[BestOverall],
                       "\n---------------------------------------------------\n"))

  }
  #### Return FunctionOutput here ##
  FunctionOutput = list(message = message, FunctionOutput = FunctionOutput)
  class(FunctionOutput) <- c("ClusterVARSummary", class(FunctionOutput))

  return(FunctionOutput)


} # eoF
