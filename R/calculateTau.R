calculateTau <- function(tau, FZY, N) {
    tau <- colSums(FZY) / N
    invisible(tau)
}