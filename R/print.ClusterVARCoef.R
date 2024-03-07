print.ClusterVARCoef <- function(x, ...) {
    #### Return FunctionOutput here ##
    FunctionOutput = list(Lags = x$Lags,
                          VAR_coefficients = x$VAR_coefficients,
                          Exogenous_coefficients = x$Exogenous_coefficients,
                          Sigma = x$Sigma,
                          Proportions = x$Proportions,
                          PredictableTimepoints = x$PredictableTimepoints,
                          Converged = x$Converged)
    print(FunctionOutput)
} # eoF
