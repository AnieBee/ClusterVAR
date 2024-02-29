
coef.ClusterVAR <- function(object, ...) {

  # Collect from ...
  args <- list(...)
  if(is.null("args$Model")) stop("Specify models via the Models argument.")
  Lags <- args$Model

  # Lags: a numeric integer or vector of length equal to the number of clusters in the model of interest (do not need to be subsequent).
  # Specifies the lag order of the model for which model coefficients should be retrieved
  # Returns only the important coefficients of the best start for a given model (i.e., given the lag order and number of clusters)

  #get best model for this lag order and this number of clusters
  nClusters = length(Lags)
  # -------------Check all input is as expected---------------
  if(!(nClusters %in% object$Call$Clusters)){
    stop("The length of 'Models', which specifies the number of clusters, is not contained in the Cluster Sequence of your fitted LCVAR Model Object.")
  }
  if(!all(Lags %in% object$Call$Lags)){
    stop("Models', which specifies the lag order, is not contained in the Lag Sequence of your fitted LCVAR Model Object.")
  }
  # -------------

  ### Find chunk of output that contains all starts for the correct lag order and number of clusters ###
  Lags = Lags[order(Lags)]
  LagListClust = calculateLagList(K = nClusters, HighestLag = max(object$Call$Lags), LowestLag = min(object$Call$Lags))
  Lagscounter  = which(apply(LagListClust, 1, function(row) all(row == Lags)))
  clusterCounter = which(object$Call$Clusters == nClusters)
  GivenOutput = object$All_Models[[clusterCounter]][[Lagscounter]]

  ## Start loop to find start with best loglikelihood value
  NumberStarts = (object$Call$Rand + as.numeric(object$Call$Rational) + as.numeric(!is.null(object$Call$Initialization)) + 1)
  FitAllStarts = array(NA, dim = NumberStarts) # to store fit of output

  for(StartCounter in 1:NumberStarts){
    FitAllStarts[StartCounter] = GivenOutput[[StartCounter]]$last.loglik
  }
  best_model_index = which.max(FitAllStarts)

  FunctionOutput = GivenOutput[[best_model_index]]
  class(FunctionOutput) <- c("ClusterVARCoef", class(FunctionOutput))

  return(FunctionOutput)

} # eoF
