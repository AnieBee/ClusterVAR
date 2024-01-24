

print.ClusterVAR <- function(x) {

  cat("----- Fitted LCVAR Model Object -----")
  # cat("\n")
  cat(paste0("\nCluster Sequence: ",  paste0(x$Call$Clusters, collapse=" ")))
  cat(paste0("\nLag Sequence: ",  paste0(x$Call$Lags, collapse=" ")))
  cat("\nRuntime: ",  x$Runtime)

} # eoF
