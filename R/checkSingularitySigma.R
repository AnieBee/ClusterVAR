checkSingularitySigma <- function(nDepVar, K, Sigma) {
    iterationReset <- FALSE

    # Check for singularity using a vectorized approach
    singular_components <- sapply(1:K, function(j) det(Sigma[, , j]) < 1.0e-200)

    # Identify the singular components
    singular_indices <- which(singular_components)

    # Reset the diagonals for singular components
    if(any(singular_components)) Sigma[, , singular_indices] <- Sigma[, , singular_indices] + 0.01

    # Update iterationReset flag
    iterationReset <- any(singular_components)

    invisible(list(Sigma = Sigma, iterationReset = iterationReset))
}
