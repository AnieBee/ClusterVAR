#
# test_that("mlVAR Perm Test check: P-values", {
#
#   # Reproducibility
#   set.seed(1)
#
#   # Call Permutation test on exampe data
#   out <- mlVAR_GC(data = ExampleData,
#                   vars = c("V1", "V2", "V3"),
#                   idvar = "id",
#                   groups = "group",
#                   nCores = 1,
#                   nP = 5)
#
#   # Check whether group difference is the same
#   expect_equal(out$Pval$Lagged_fixed[1,2],  0)
#
# })
#
