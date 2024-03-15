# Somtimes but not always faster than non-parallelized B
# Not sure if the list that is returned is always in the same order in terms of 1:K, B[ , , j] has to return the B for the correct cluster
calculateB <- function(Covariates, K, nDepVar, A, Sigma, N, NewPredictableObs, X, Y, Lags, FZY, qqq, B)
{
 ############NEW
  print("Loading ClusterVAR package...")
  devtools::load_all("/home/anja/Desktop/ClusterVARPackage/ClusterVAR")
  print("ClusterVAR package loaded successfully.")
  ###########

    #if(Covariates == "equal-within-clusters"){# B(0): equal-within-clusters
        cluster <- 6
        cl <- makeCluster(cluster, outfile="")
        registerDoParallel(cl)

        B_List <- foreach(j = 1:K,
                .packages = c("MASS")) %dopar%
        {

            Bnum = 0
            Bdenom = 0
            Delta = cbind(diag(nDepVar), (-1) * A[ , 1:(nDepVar * Lags[j]), j])
            NewCovaMat = t(Delta) %*% MASS::ginv(Sigma[, , j]) %*% Delta
            for(i in 1:N)
            {
                Bn = 0
                Bd = 0
                for(trunner in c(NewPredictableObs[[ Lags[j] ]][[i]]) )
                {
                    XtildaKron <- kronecker(t(X[ , (trunner):(trunner - Lags[j]), drop = FALSE]),
                                            diag(nDepVar)) # drop = FALSE is needed here in case of q = 1
                    Bd = Bd + ( t(XtildaKron) %*% NewCovaMat %*% XtildaKron )
                    Bn = Bn + ( t(XtildaKron) %*% NewCovaMat %*%
                                    as.vector(Y[ , (trunner):(trunner - Lags[j]), drop = FALSE]) )
                }
                Bnum = Bnum + (FZY[ i, j] * Bn)
                Bdenom = Bdenom + (FZY[ i, j] * Bd)
            }
            BasVec = MASS::ginv(Bdenom) %*% Bnum
            return(matrix(BasVec, nrow = nDepVar, ncol = qqq, byrow = FALSE))
        }
    #}
    stopCluster(cl)
    for(j in 1:K) B[, , j] = B_List[[j]]

    invisible(B)

}
