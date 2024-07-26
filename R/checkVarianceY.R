checkVarianceY <- function(Y, NewPredictableObs, N, HighestLag, pers)
{
    PrintInfoWarning = FALSE
    for(i in 1:N){
        LogicalZeroVariance = apply(Y[ , NewPredictableObs[[HighestLag]][[i]], drop = FALSE], 1, var) < .00001
        if(any(LogicalZeroVariance)){
            warning(paste("For the highest number of 'Lags', individual", pers[i],"has virtually no variance in predictable observations for the variable:", names(LogicalZeroVariance)[which(LogicalZeroVariance)], "\n"))
            PrintInfoWarning = TRUE
        } 
    }
    
    if(PrintInfoWarning) warning(paste("Including individuals with no variance in the endogenous variables can lead to pathological solutions, particularly when there are many such cases. Please consider removing these individual(s) or variable(s). See the paper for details."))
}