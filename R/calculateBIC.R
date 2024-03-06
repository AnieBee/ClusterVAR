calculateBIC = function(nPara, clTimepoints, last.lik)
{
    invisible( (-2 * last.lik) + (nPara * log(sum(clTimepoints))) )
}