
summary.ClusterVAR <- function(CVARresult, show = "Best-per-number-of-clusters", TS_criterion = "SC", global_criterion = "BIC",
                               Number_of_Clusters = NULL) {
    #show = c("Best-per-number-of-clusters", "Best-overall", "Given-a-number-of-clusters")
    #TS_criterion = c("SC", "HQ")
    #Global_criterion = c("BIC", "ICL")
    # Number_of_Clusters is only used if show == "Given-a-number-of-clusters"
    
    # -------------Check all input is as expected---------------
    if (!(show %in% c("Best-per-number-of-clusters", "Best-overall", "Given-a-number-of-clusters"))) {
        stop("Invalid value for 'show'. Please choose from: ", paste(c("Best-per-number-of-clusters", "Best-overall", "Given-a-number-of-clusters"), collapse = ", "))
    }
    if (!(TS_criterion %in% c("SC", "HQ"))) {
        stop("Invalid value for 'TS_criterion'. Please choose from: ", paste(c("SC", "HQ"), collapse = ", "))
    }
    if (!(global_criterion %in% c("BIC", "ICL"))) {
        stop("Invalid value for 'global_criterion'. Please choose from: ", paste(c("BIC", "ICL"), collapse = ", "))
    }
    if (show == "Given-a-number-of-clusters" && is.null(Number_of_Clusters)) {
        stop("If 'show' is 'Given-a-number-of-clusters', you must specify a value for 'Number_of_Clusters'.")
    }
    if(show == "Given-a-number-of-clusters" && !(Number_of_Clusters %in% CVARresult$Call$Clusters)){
        stop("The value you specified for 'Number_of_Clusters' is not contained in the Cluster Sequence of your fitted LCVAR Model Object.")
    }
    # -------------
    
    l_LagsPerCluster <- list()
    for(K in CVARresult$Call$Clusters)  l_LagsPerCluster[[K]] <- calculateLagList(K = K, HighestLag = max(CVARresult$Call$Lags), LowestLag = min(CVARresult$Call$Lags))
    NumberStarts = (CVARresult$Call$Rand + as.numeric(CVARresult$Call$Rational) + as.numeric(!is.null(CVARresult$Call$Initialization)) + 1) 
    # last element is for PreviousSol which is TRUE as a default
    
    # -------------Create one of the three different types of FunctionOutput---------------
    
    if(show == "Given-a-number-of-clusters"){
        
        clusterCounter = which(CVARresult$Call$Clusters == Number_of_Clusters)
        GivenClusterOutput = CVARresult$All_Solutions_All_Clusters_Lags_Starts[[clusterCounter]]
        LagCombinations = nrow(l_LagsPerCluster[[Number_of_Clusters]])
        
        # Fit[Lags, Start]
        FitAllLags = array(NA, dim = c(LagCombinations, NumberStarts)) # to store fit of output
        FunctionOutput = data.frame(matrix(NA, nrow = LagCombinations, ncol = 6), 
                                    row.names = apply(l_LagsPerCluster[[Number_of_Clusters]], 1, function(x) paste(c("Lags", x), collapse = " ")))
        colnames(FunctionOutput) = c("log-likelihood", "parameters", "SC", "HQ", "BIC", "ICL")
       
        
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
        }
        cat(paste0(c("--------------------------------------------------------------------------------------------------------", "\n",
                     "All lags for number of clusters =", Number_of_Clusters,
                     "\n",
                     "--------------------------------------------------------------------------------------------------------\n")))
    }else{
        # The below calculation calculates the solution for (show == "Best-per-number-of-clusters") but all calculation steps are also needed if (show == "Best-overall") 
            # For each number of of clusters, find the single best-fitting time-series model for each cluster number
           
            FunctionOutput = data.frame(matrix(NA, nrow = length(CVARresult$Call$Clusters), ncol = 5), 
                                        row.names = apply(as.matrix(CVARresult$Call$Clusters), 1, function(x) paste(c(x, "Clusters"), collapse = " ")))
            colnames(FunctionOutput) = c(paste(c("Lags selected by", TS_criterion), collapse = " "), "log-likelihood", "parameters", "BIC", "ICL")
            
            ClustCount = 1
            for(K in CVARresult$Call$Clusters){
                
                LagCombinations = nrow(l_LagsPerCluster[[K]])
                FitAllLags = array(NA, dim = c(LagCombinations, NumberStarts)) # to store fit of output
                
                for(LagCounter in 1:LagCombinations) {
                    for(StartCounter in 1:NumberStarts){
                        FitAllLags[LagCounter, StartCounter] = switch(TS_criterion,
                                                                      "SC" = CVARresult$All_Solutions_All_Clusters_Lags_Starts[[ClustCount]][[LagCounter]][[StartCounter]]$SC,
                                                                      "HQ" = CVARresult$All_Solutions_All_Clusters_Lags_Starts[[ClustCount]][[LagCounter]][[StartCounter]]$HQ)
                    }
                }
                BestRunOneK = arrayInd(which.min(FitAllLags), dim(FitAllLags))
                FunctionOutput[ClustCount, 1] = paste(CVARresult$All_Solutions_All_Clusters_Lags_Starts[[ClustCount]][[BestRunOneK[1]]][[BestRunOneK[2]]]$Lags, collapse = " ")
                FunctionOutput[ClustCount, "log-likelihood"] = CVARresult$All_Solutions_All_Clusters_Lags_Starts[[ClustCount]][[BestRunOneK[1]]][[BestRunOneK[2]]]$last.loglik
                FunctionOutput[ClustCount, "parameters"] = CVARresult$All_Solutions_All_Clusters_Lags_Starts[[ClustCount]][[BestRunOneK[1]]][[BestRunOneK[2]]]$nParameters
                FunctionOutput[ClustCount, "BIC"] = CVARresult$All_Solutions_All_Clusters_Lags_Starts[[ClustCount]][[BestRunOneK[1]]][[BestRunOneK[2]]]$BIC
                FunctionOutput[ClustCount, "ICL"] = CVARresult$All_Solutions_All_Clusters_Lags_Starts[[ClustCount]][[BestRunOneK[1]]][[BestRunOneK[2]]]$ICL
                
                ClustCount = ClustCount + 1
            }
            
            
            if(show == "Best-per-number-of-clusters"){
                cat(paste0(c("--------------------------------------------------------------------------------------------------------", "\n",
                             "The best lags for any number of clusters as selected by the", TS_criterion,
                         "\n",
                         "--------------------------------------------------------------------------------------------------------\n")))
                # FunctionOutput stays as calculated above
            }
            if (show == "Best-overall"){
                BestOverall = switch(global_criterion,
                                    "BIC" = which.min(FunctionOutput$BIC),
                                    "ICL" = which.min(FunctionOutput$ICL))
                cat(paste0(c("--------------------------------------------------------------------------------------------------------", "\n",
                             "The best lags for any number of clusters as selected by the", TS_criterion,
                             "\n",
                             "The best number of clusters as selected by the", global_criterion,
                             "\n",
                             "-------------------------------------------------------------------------------------------------------- \n")))
                FunctionOutput = FunctionOutput[BestOverall, ]
            }

    }
    #### Return FunctionOutput here ##
    return(FunctionOutput)
} # eoF
