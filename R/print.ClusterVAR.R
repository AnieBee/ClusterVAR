

print.ClusterVAR <- function(x) {

  cat("----- Fitted LCVAR Model Object -----")
  # cat("\n")
  cat(paste0("\nClusters Sequence: ",  paste0(x$Call$Clusters, collapse=" ")))
  cat(paste0("\nLags Sequence: ",  paste0(x$Call$Lags, collapse=" ")))
  cat("\nRuntime: ",  x$Runtime, " mins")

} # eoF
