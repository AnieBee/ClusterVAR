calculateICL = function(BIC, K, N, FZY, Classification)
{
    ICLPenaltyTerm = 0
    for(j in 1:K)
    {
        for (i in 1:N)
        {
            ICLPenaltyTerm = ICLPenaltyTerm + ( as.numeric(j == Classification[i]) * ifelse(FZY[ i, j] == 0, log(1e-260), log(FZY[ i, j])) )
            # in case the posterior prob of cluster membership (FZY) is zero >> replace with 1-260 to avoid getting -Inf as a result of log(0)
        }
    }
    
    invisible( (BIC - (2 * ICLPenaltyTerm)) )
}


