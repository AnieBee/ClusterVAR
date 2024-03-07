calculateBandWZero <- function(Covariates, K, memb, Y, X , PersStart, PersEnd, BZero, WZero)
{
    
    if(Covariates == "equal-within-clusters"){ # Same within clusters, dim(B)[3] = K
        for(j in 1:K)
        {
            BZeronum <- 0 
            BZerodenom <- 0 
            for(i in c(which(memb[j, ] == 1)))
            {  
                XPers = X[ , (PersStart[i]):(PersEnd[i]), drop = FALSE]
                BZeronum <- BZeronum + Y[ , PersStart[i]:PersEnd[i], drop = FALSE] %*% t(XPers)
                BZerodenom <- BZerodenom + XPers %*% t(XPers)
            }
            BZero[ , , j] <- BZeronum%*%MASS::ginv(BZerodenom)
            for(i in c(which(memb[j, ] == 1)))
            {  
                WZero[ , (PersStart[i]):(PersEnd[i]), 1] <- Y[ , (PersStart[i]):(PersEnd[i]), drop = FALSE] - (BZero[ , , j] %*% 
                                                            X[ , (PersStart[i]):(PersEnd[i]), drop = FALSE])
            }
        }
    }#else{# Anja if B is constraint, Wk is same for all clusters and B is same for all clusters  dim(B)[3] = 1
    #     if(Covariates == "equal-across-clusters")
    #     { 
    #         BZero[ , , 1] <- (Y %*% t(X)) %*% ginv(X %*% t(X))
    #         WZero[ , , 1] <- Y - BZero[ , , 1] %*% X
    #     } # else
    #     # { 
    #     # }
    # }
    
    invisible(WZero)
    
}

