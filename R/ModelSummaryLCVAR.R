ModelSummaryLCVAR <- function(LCVARclustResult, Output = c("LogLikelihood", "IC"))
{    
    # Shall we include a check here that the LCVARclustResult input really is output from the LCVARclust function?
    
    if(is.null(Output)){Output = "LogLikelihood"}
    if(Output != "LogLikelihood" & Output != "IC"){Output = "LogLikelihood"}
  
    ModelCall = LCVARclustResult$Call
    Summary = vector("list", length(ModelCall$Clusters))
    #Create a name for all the list elements of summary (for each cluster)
    listnames = NULL
    OutputName = ifelse(Output == "LogLikelihood", "Log likelihood",
                        paste(c("Information criteria", ModelCall$ICType), collapse = " "))
    
    ClusterCounter = 0
    for (K in ModelCall$Clusters)
    {
        ClusterCounter = ClusterCounter + 1
        
        listnames = c(listnames,
                      paste(c(OutputName, "for all models where the number of clusters is:", ModelCall$Clusters[ClusterCounter]), collapse = " "))
        
        LagsList = calculateLagList(K = K, HighestLag = ModelCall$HighestLag, LowestLag = ModelCall$LowestLag)
        if (K == 1){ LagsList = t(LagsList) } # dimension changes when number of clusters is equal to 1 >> In that case, transpose to get same dimension again
        LagCombinations = dim(LagsList)[1]
        
        PreviousSol = TRUE 
        AllModelsThisCluster = array(NA, dim = c(LagCombinations, ModelCall$Rand + as.numeric(ModelCall$Rational) + 
                                       as.numeric(!is.null(ModelCall$Initialization))
                                        + as.numeric(PreviousSol))) # to store fit of output
        
        ## Name Columns and rows for the output
        colnames(AllModelsThisCluster) = c(paste(c(paste("Random Start Nr. ", as.character(1:ModelCall$Rand), ":" )), sep = ""),
                      ifelse(ModelCall$Rational, "Rational Start:", NULL),
                      switch(!is.null(ModelCall$Initialization), "FALSE" = "Initialization Start:"), # default is NULL
                      ifelse(PreviousSol, "Previous Start:", NULL))
        rowname = NULL
        for (i in 1:nrow(AllModelsThisCluster)){
            rowname = c(rowname, paste(c(LagsList[i, ], "Lags", ":"), collapse = " "))
        }
        rownames(AllModelsThisCluster) = rowname
        
        if(Output == "LogLikelihood"){
            for (LagCounter in 1:LagCombinations){
                for(StartCounter in 1:(dim(AllModelsThisCluster)[2])){
                    AllModelsThisCluster[LagCounter,  StartCounter] =  LCVARclustResult$AllSolutions[[ClusterCounter]][[LagCounter]][[StartCounter]]$last.loglik
                }
            }
        }
        if(Output == "IC"){
            for (LagCounter in 1:LagCombinations){
                for(StartCounter in 1:(dim(AllModelsThisCluster)[2])){
                    AllModelsThisCluster[LagCounter,  StartCounter] =  LCVARclustResult$AllSolutions[[ClusterCounter]][[LagCounter]][[StartCounter]]$IC
                }
            }
        }
        
        Summary[[ClusterCounter]] = t(AllModelsThisCluster)

    }
    
    names(Summary) = listnames
    
    return(Summary)
}