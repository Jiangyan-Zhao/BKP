# This test file validates the functionality of the `parameter` S3 method for
# DKP objects. It ensures that the method correctly extracts the posterior
# Dirichlet parameters and the optimized kernel hyperparameters from a fitted model.

test_that("parameter returns correct parameters for DKP models", {
  # Set a seed for reproducibility
  set.seed(42)

  # Generate synthetic data for testing
  d <- 2
  n <- 50
  q <- 3 # number of classes
  X <- matrix(runif(n * d), n, d)
  m <- rep(10, n)
  Y <- t(sapply(1:n, function(i) {
    p <- c(0.2, 0.3, 0.5)
    rmultinom(1, size = m[i], prob = p)
  }))

  # -------------------------------------------------------------------------
  # Test Case 1: Noninformative Prior
  # -------------------------------------------------------------------------

  # Fit the DKP model
  model_noninformative <- fit_DKP(X, Y, prior = "noninformative")

  # Extract parameters using the S3 method
  params_noninformative <- parameter(model_noninformative)

  # Verify that the extracted parameters match the values stored in the model object
  expect_equal(params_noninformative$theta, model_noninformative$theta_opt, tolerance = 1e-6)
  expect_equal(params_noninformative$alpha_n, model_noninformative$alpha_n, tolerance = 1e-6)

  # -------------------------------------------------------------------------
  # Test Case 2: Fixed Prior
  # -------------------------------------------------------------------------

  r0 <- 10
  p0 <- c(0.25, 0.5, 0.25)

  # Fit the DKP model
  model_fixed <- fit_DKP(X, Y, prior = "fixed", r0 = r0, p0 = p0)

  # Extract parameters using the S3 method
  params_fixed <- parameter(model_fixed)

  # Verify that the extracted parameters match the values stored in the model object
  expect_equal(params_fixed$theta, model_fixed$theta_opt, tolerance = 1e-6)
  expect_equal(params_fixed$alpha_n, model_fixed$alpha_n, tolerance = 1e-6)

  # -------------------------------------------------------------------------
  # Test Case 3: Adaptive Prior
  # -------------------------------------------------------------------------

  r0 <- 10

  # Fit the DKP model
  model_adaptive <- fit_DKP(X, Y, prior = "adaptive", r0 = r0)

  # Extract parameters using the S3 method
  params_adaptive <- parameter(model_adaptive)

  # Verify that the extracted parameters match the values stored in the model object
  expect_equal(params_adaptive$theta, model_adaptive$theta_opt, tolerance = 1e-6)
  expect_equal(params_adaptive$alpha_n, model_adaptive$alpha_n, tolerance = 1e-6)
})
