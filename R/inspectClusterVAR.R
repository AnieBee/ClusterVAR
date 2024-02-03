inspectClusterVAR <- function(CVARresult, Lags) {
    #Lags: a numeric integer or vector of length equal to the number of clusters in the model of interest (do not need to be subsequent).
    #Specifies the lag order of the model for which model coefficients should be retrieved
    # Returns all info of the best start for a given model (i.e., given the lag order and number of clusters)
    nClusters = length(Lags)
    # -------------Check all input is as expected---------------
    if(!(nClusters %in% CVARresult$Call$Clusters)){
        stop("The length of 'Lags', which specifies the number of clusters, is not contained in the Cluster Sequence of your fitted LCVAR Model Object.")
    }
    if(!all(Lags %in% CVARresult$Call$Lags)){
        stop("Lags', which specifies the lag order, is not contained in the Lag Sequence of your fitted LCVAR Model Object.")
    }
    # -------------
    
    ### Find chunk of output that contains all starts for the correct lag order and number of clusters ###
    Lags = Lags[order(Lags)]
    LagListClust = calculateLagList(K = nClusters, HighestLag = max(CVARresult$Call$Lags), LowestLag = min(CVARresult$Call$Lags))
    Lagscounter  = which(apply(LagListClust, 1, function(row) all(row == Lags)))
    clusterCounter = which(CVARresult$Call$Clusters == nClusters)
    GivenOutput = CVARresult$All_Solutions_All_Clusters_Lags_Starts[[clusterCounter]][[Lagscounter]]
    
    ## Start loop to find start with best loglikelihood value
    NumberStarts = (CVARresult$Call$Rand + as.numeric(CVARresult$Call$Rational) + as.numeric(!is.null(CVARresult$Call$Initialization)) + 1) 
    FitAllStarts = array(NA, dim = NumberStarts) # to store fit of output
    for(StartCounter in 1:NumberStarts){
        FitAllStarts[StartCounter] = GivenOutput[[StartCounter]]$last.loglik
    }
    best_model_index = which.max(FitAllStarts)
    
    #### Return FunctionOutput here ##
    return(GivenOutput[[best_model_index]])
    
} # eoF