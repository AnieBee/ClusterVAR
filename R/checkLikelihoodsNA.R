checkLikelihoodsNA <- function(FYZ, EMiteration)
{
    iterationReset = FALSE
    stopifnot(sum(FYZ, na.rm = TRUE) != 0) # check for underflow, if all posteriors are zero
    if(any(is.na(FYZ)))
    { # Check not a single posterior is NA
        if(all(is.na(FYZ)))
        {
            stop("All likelihoods are NA in EM-iteration", EMiteration, "\n")
        } 
        #cat("\n Warning: Some likelihoods NA in EM-iteration", EMiteration, ", likelihoods reset \n")
        NAIndexFYZ = which(is.na(FYZ))
        FYZ[NAIndexFYZ] = mean(FYZ, na.rm=TRUE) # Set those p(Y|z_{ik}) that are NA for some person to the overall mean of all P(Y|Z_{ik})
        iterationReset = TRUE
    }
    
    invisible(list(FYZ = FYZ, iterationReset = iterationReset))
}