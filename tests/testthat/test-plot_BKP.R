# This test file validates the functionality of the `plot` S3 method for
# BKP objects. It ensures that the function generates plots without
# errors for different input dimensions and parameter settings.

test_that("plot.BKP generates plots without errors", {
  # Set a seed for reproducibility
  set.seed(123)

  # -------------------------------------------------------------------------
  # Test Case 1: 1D Input
  # -------------------------------------------------------------------------

  d <- 1
  n <- 50
  X <- matrix(runif(n * d), n, d)
  y <- rbinom(n, size = 10, prob = 0.5)
  m <- rep(10, n)

  # Fit a 1D BKP model
  model_1d <- fit_BKP(X, y, m, prior = "noninformative")

  # Test that plot() runs without errors for a 1D model
  expect_no_error(plot(model_1d))

  # -------------------------------------------------------------------------
  # Test Case 2: 2D Input
  # -------------------------------------------------------------------------

  d <- 2
  n <- 50
  X <- matrix(runif(n * d), n, d)
  y <- rbinom(n, size = 10, prob = 0.5)
  m <- rep(10, n)

  # Fit a 2D BKP model
  model_2d <- fit_BKP(X, y, m, prior = "fixed", r0 = 10, p0 = 0.5)

  # Test with default arguments
  expect_no_error(plot(model_2d, n_grid = 30))

  # Test with only_mean = TRUE
  expect_no_error(plot(model_2d, only_mean = TRUE, n_grid = 30))

  # Test with a smaller n_grid
  expect_no_error(plot(model_2d, n_grid = 30))
  # -------------------------------------------------------------------------
  # Test Case 3: 3D Input
  # -------------------------------------------------------------------------

  d <- 3
  n <- 50
  X <- matrix(runif(n * d), n, d)
  y <- rbinom(n, size = 10, prob = 0.5)
  m <- rep(10, n)

  # Fit a 2D BKP model
  model_2d <- fit_BKP(X, y, m, prior = "fixed", r0 = 10, p0 = 0.5)

  # Test with default arguments
  expect_no_error(plot(model_2d, dims = c(1,3), n_grid = 10))

  # Test with only_mean = TRUE
  expect_no_error(plot(model_2d, only_mean = TRUE, dims = c(1,3), n_grid = 10))
})
