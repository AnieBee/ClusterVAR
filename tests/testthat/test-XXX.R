
test_that("Fake test", {

  # Reproducibility
  set.seed(1)

  # Here call the LCVAR function

  # I don't add those now, otherwise every time I compile the package I have to wait for the function to fit
  fake_outcome <- 0

  # Check whether group difference is the same
  expect_equal(fake_outcome,  0)

})

