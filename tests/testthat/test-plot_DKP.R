# This test file validates the functionality of the `plot` S3 method for
# DKP objects. It ensures that the function generates plots without
# errors for 1D and 2D input dimensions.

test_that("plot.DKP generates plots without errors for 1D and 2D inputs", {
  # Set a seed for reproducibility to ensure consistent data
  set.seed(123)

  # -------------------------------------------------------------------------
  # Test Case 1: 1D Input
  # -------------------------------------------------------------------------

  # Define a spatially-varying probability function for 1D data.
  # This function must return a vector of length 'q' for each input point.
  true_pi_fun_1d <- function(x) {
    # Define probabilities for each class
    p1 <- exp(sin(2 * pi * x))
    p2 <- exp(sin(4 * pi * x))
    p3 <- exp(cos(2 * pi * x))

    # Normalize to ensure probabilities sum to 1
    total <- p1 + p2 + p3
    return(c(p1 / total, p2 / total, p3 / total))
  }

  d <- 1
  n <- 50
  q <- 3 # number of classes
  X <- matrix(runif(n * d), n, d)
  m <- rep(10, n)

  Y <- t(sapply(1:n, function(i) {
    p <- true_pi_fun_1d(X[i])
    rmultinom(1, size = m[i], prob = p)
  }))

  # Fit a 1D DKP model
  model_1d <- fit_DKP(X, Y, prior = "noninformative")

  # Test that plot() runs without errors for a 1D model
  expect_no_error(plot(model_1d))

  # -------------------------------------------------------------------------
  # Test Case 2: 2D Input
  # -------------------------------------------------------------------------

  # Define a spatially-varying probability function for 2D data
  true_pi_fun_2d <- function(X) {
    p1 <- exp(sin(2 * pi * X[,1]) + cos(2 * pi * X[,2]))
    p2 <- exp(sin(2 * pi * X[,2]) - cos(2 * pi * X[,1]))
    p3 <- exp(sin(4 * pi * X[,1]))

    # Normalize to ensure probabilities sum to 1
    total <- p1 + p2 + p3
    return(cbind(p1/total, p2/total, p3/total))
  }

  d <- 2
  n <- 50
  q <- 3
  X <- matrix(runif(n * d), n, d)
  m <- rep(10, n)

  Y <- t(sapply(1:n, function(i) {
    p <- true_pi_fun_2d(X[i, , drop = FALSE])
    rmultinom(1, size = m[i], prob = p)
  }))

  # Fit a 2D DKP model
  model_2d <- fit_DKP(X, Y, prior = "fixed", r0 = 10, p0 = c(1/3, 1/3, 1/3))

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
  Y <- t(rmultinom(n, size = 10, prob = c(1/3,1/3,1/3)))
  # Fit a 2D DKP model
  model_2d <- fit_DKP(X, Y, prior = "fixed", r0 = 10, p0 = c(1/3, 1/3, 1/3))

  # Test with default arguments
  expect_no_error(plot(model_2d, dims=c(1,3), n_grid = 10))

  # Test with only_mean = TRUE
  expect_no_error(plot(model_2d, only_mean = TRUE, dims=c(1,3), n_grid = 10))
})
