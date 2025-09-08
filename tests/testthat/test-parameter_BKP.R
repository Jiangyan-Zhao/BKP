# This test file validates the functionality of the `parameter` S3 method for
# BKP objects. It ensures that the method correctly extracts the posterior
# Beta parameters and the optimized kernel hyperparameters from a fitted model.

test_that("parameter returns correct parameters for BKP models", {
  # Set a seed for reproducibility
  set.seed(42)

  # Generate synthetic data for testing
  d <- 2
  n <- 50
  X <- matrix(runif(n * d), n, d)
  y <- rbinom(n, size = 10, prob = 0.5)
  m <- rep(10, n)

  # -------------------------------------------------------------------------
  # Test Case 1: Noninformative Prior
  # -------------------------------------------------------------------------

  # Fit the BKP model
  model_noninformative <- fit_BKP(X, y, m, prior = "noninformative")

  # Extract parameters using the S3 method
  params_noninformative <- parameter(model_noninformative)

  # Verify that the extracted parameters match the values stored in the model object
  expect_equal(params_noninformative$alpha_n, model_noninformative$alpha_n, tolerance = 1e-6)
  expect_equal(params_noninformative$beta_n, model_noninformative$beta_n, tolerance = 1e-6)

  # -------------------------------------------------------------------------
  # Test Case 2: Fixed Prior
  # -------------------------------------------------------------------------

  r0 <- 10
  p0 <- 0.6

  # Fit the BKP model
  model_fixed <- fit_BKP(X, y, m, prior = "fixed", r0 = r0, p0 = p0)

  # Extract parameters using the S3 method
  params_fixed <- parameter(model_fixed)

  # Verify that the extracted parameters match the values stored in the model object
  expect_equal(params_fixed$alpha_n, model_fixed$alpha_n, tolerance = 1e-6)
  expect_equal(params_fixed$beta_n, model_fixed$beta_n, tolerance = 1e-6)

  # -------------------------------------------------------------------------
  # Test Case 3: Adaptive Prior
  # -------------------------------------------------------------------------

  r0 <- 10

  # Fit the BKP model
  model_adaptive <- fit_BKP(X, y, m, prior = "adaptive", r0 = r0)

  # Extract parameters using the S3 method
  params_adaptive <- parameter(model_adaptive)

  # Verify that the extracted parameters match the values stored in the model object
  expect_equal(params_adaptive$alpha_n, model_adaptive$alpha_n, tolerance = 1e-6)
  expect_equal(params_adaptive$beta_n, model_adaptive$beta_n, tolerance = 1e-6)
})
