coef.ClusterVAR <- function(LCVARresult, Lags) {

    # Lags: a numeric integer or vector of length equal to the number of clusters in the model of interest (do not need to be subsequent).
    # Specifies the lag order of the model for which model coefficients should be retrieved
    # Returns only the important coefficients of the best start for a given model (i.e., given the lag order and number of clusters)

    #get best model for this lag order and this number of clusters
    nClusters = length(Lags)
    # -------------Check all input is as expected---------------
    if(!(nClusters %in% LCVARresult$Call$Clusters)){
        stop("The length of 'Lags', which specifies the number of clusters, is not contained in the Cluster Sequence of your fitted LCVAR Model Object.")
    }
    if(!all(Lags %in% LCVARresult$Call$Lags)){
        stop("Lags', which specifies the lag order, is not contained in the Lag Sequence of your fitted LCVAR Model Object.")
    }
    # -------------
    
    ### Find chunk of output that contains all starts for the correct lag order and number of clusters ###
    Lags = Lags[order(Lags)]
    LagListClust = calculateLagList(K = nClusters, HighestLag = max(LCVARresult$Call$Lags), LowestLag = min(LCVARresult$Call$Lags))
    Lagscounter  = which(apply(LagListClust, 1, function(row) all(row == Lags)))
    clusterCounter = which(LCVARresult$Call$Clusters == nClusters)
    GivenOutput = LCVARresult$All_Solutions_All_Clusters_Lags_Starts[[clusterCounter]][[Lagscounter]]
    
    ## Start loop to find start with best loglikelihood value
    NumberStarts = (LCVARresult$Call$Rand + as.numeric(LCVARresult$Call$Rational) + as.numeric(!is.null(LCVARresult$Call$Initialization)) + 1)
    FitAllStarts = array(NA, dim = NumberStarts) # to store fit of output
    for(StartCounter in 1:NumberStarts){
        FitAllStarts[StartCounter] = GivenOutput[[StartCounter]]$last.loglik
    }
    best_model_index = which.max(FitAllStarts)
    
    FunctionOutput = GivenOutput[[best_model_index]]
    class(FunctionOutput) <- c("ClusterVARCoef", class(FunctionOutput))
    
    return(FunctionOutput)

} # eoF
