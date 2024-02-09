coef.ClusterVAR <- function(CVARresult, Lags) {

    # Lags: a numeric integer or vector of length equal to the number of clusters in the model of interest (do not need to be subsequent).
    # Specifies the lag order of the model for which model coefficients should be retrieved
    # Returns only the important coefficients of the best start for a given model (i.e., given the lag order and number of clusters)

    #get best model for this lag order and this number of clusters
    GivenOutput = inspectClusterVAR(CVARresult, Lags = Lags)

    #### Return FunctionOutput here ##
    FunctionOutput = list(Lags = GivenOutput$Lags,
                          VAR_coefficients = GivenOutput$VAR_coefficients,
                          Exogenous_coefficients = GivenOutput$Exogenous_coefficients,
                          Sigma = GivenOutput$Sigma,
                          Proportions = GivenOutput$Proportions)
    return(FunctionOutput)

} # eoF
