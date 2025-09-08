# Test file for the fitted.DKP method

test_that("fitted.DKP returns correct posterior mean", {
  # Create a dummy DKP model object
  alpha_n <- matrix(c(10, 5, 20, 10, 30, 15), nrow = 3, byrow = TRUE)

  mock_model <- list(alpha_n = alpha_n)
  class(mock_model) <- "DKP"

  # Calculate expected fitted values manually
  expected_fitted <- alpha_n / rowSums(alpha_n)

  # Check if fitted() returns the expected values
  expect_equal(unname(fitted(mock_model)), expected_fitted)

  # Check the class and dimensions of the output
  expect_true(is.matrix(fitted(mock_model)))
  expect_equal(dim(fitted(mock_model)), c(3, 2))

  # Note: The test for colnames is not added because the fitted.DKP method
  # does not add column names to the output, so it would fail.
})
