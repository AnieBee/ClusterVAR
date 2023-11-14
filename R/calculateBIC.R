calculateBIC = function(nPara, Lags, K, N, FZY, Tni, last.lik)
{
    clTimepoints = as.vector(rep(0, K))
    for(j in 1:K)
    {
        for (i in 1:N)
        {
        clTimepoints[j] = clTimepoints[j] + (FZY[ i, j] * (Tni[[Lags[j]]][i]))
        }
    }
    
    invisible( (-2 * last.lik) + (nPara * log(sum(clTimepoints))) )
}