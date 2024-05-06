
test_that("ClusterVAR", {

  # Reproducibility

  # Here call the LCVAR function
  LCVAR_outExample2 <- LCVAR(Data = ExampleData,
                             yVars = 1:4, ID = 5, Beep = 6,
                             xContinuous = 7, xFactor = 8,
                             Clusters = 2, Lags = 1, smallestClN = 3,
                             Cores = 1, RndSeed = 3, Rand = 2, Rational = TRUE,
                             SigmaIncrease = 10, it = 1, Conv = 1e-03,
                             Covariates = "equal-within-clusters")
  outcome <- coef(LCVAR_outExample2, Model = c(1, 1))
  # I don't add those now, otherwise every time I compile the package I have to wait for the function to fit
  values <- c(0.0847958, -0.1033596, 0.1068198, 0.2466590,
              -0.2176429, 0.1776731, 0.1040630, 0.1632605,
              0.2965052, -0.3046280, 0.1888623, 0.1450859,
              0.4096583, -0.4093674, 0.1878280, 0.4930010)

  # Define the row and column names
  row_names <- col_names <- c("Item1", "Item2", "Item3", "Item4")

  # Create the matrix
  matrix_data <- matrix(values, nrow = 4, ncol = 4, byrow = TRUE,
                        dimnames = list(row_names, paste0(col_names, "_t-1")))

  # Check whether group difference is the same
  expect_equal(round(outcome$VAR_coefficients[, ,1],  digits = 5),
               round(matrix_data,  digits = 5))

})

