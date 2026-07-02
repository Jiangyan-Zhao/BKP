# This test file validates the functionality of the `plot` S3 method for
# DKP objects. It ensures that the function generates plots without
# errors for 1D and 2D input dimensions.

test_that("plot.DKP generates plots without errors for 1D and 2D inputs", {
  # Set a seed for reproducibility to ensure consistent data
  set.seed(123)

  # -------------------------------------------------------------------------
  # Test Case 1: 1D Input
  # -------------------------------------------------------------------------

  set.seed(123)

  # Define true class probability function (3-class)
  true_pi_fun <- function(X) {
    p1 <- 1/(1+exp(-3*X))
    p2 <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
    return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
  }

  n <- 30
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- tgp::lhs(n = n, rect = Xbounds)
  true_pi <- true_pi_fun(X)
  m <- sample(150, n, replace = TRUE)

  # Generate multinomial responses
  Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))

  # Fit a 1D DKP model
  model_1d <- fit_DKP(X, Y, Xbounds = Xbounds, theta = 0.3, prior = "noninformative")

  # Test that plot() runs without errors for a 1D model
  expect_no_error(plot(model_1d))

  # -------------------------------------------------------------------------
  # Test Case 2: 2D Input
  # -------------------------------------------------------------------------

  set.seed(123)

  # Define latent function and transform to 3-class probabilities
  true_pi_fun <- function(X) {
    if (is.null(nrow(X))) X <- matrix(X, nrow = 1)
    m <- 8.6928; s <- 2.4269
    x1 <- 4 * X[,1] - 2
    x2 <- 4 * X[,2] - 2
    a <- 1 + (x1 + x2 + 1)^2 *
      (19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2)
    b <- 30 + (2*x1 - 3*x2)^2 *
      (18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2)
    f <- (log(a*b)- m)/s
    p1 <- pnorm(f) # Transform to probability
    p2 <- sin(pi * X[,1]) * sin(pi * X[,2])
    return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
  }

  n <- 100
  Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
  X <- tgp::lhs(n = n, rect = Xbounds)
  true_pi <- true_pi_fun(X)
  m <- sample(150, n, replace = TRUE)

  # Generate multinomial responses
  Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))

  # Fit a 2D DKP model
  model_2d <- fit_DKP(X, Y, Xbounds = Xbounds, theta = 0.3, prior = "fixed", r0 = 10, p0 = c(1/3, 1/3, 1/3))

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
  model_2d <- fit_DKP(X, Y, theta = 0.3, prior = "fixed", r0 = 10, p0 = c(1/3, 1/3, 1/3))

  # Test with default arguments
  expect_no_error(plot(model_2d, dims=c(1,3), n_grid = 10))

  # Test with only_mean = TRUE
  expect_no_error(plot(model_2d, only_mean = TRUE, dims=c(1,3), n_grid = 10))
})

test_that("plot.DKP validates input arguments and classification branches", {
  set.seed(99)

  X <- matrix(runif(60), ncol = 2)
  cl <- sample(1:3, nrow(X), replace = TRUE)
  Y <- matrix(0, nrow(X), 3)
  Y[cbind(seq_len(nrow(X)), cl)] <- 1
  model <- fit_DKP(X, Y, theta = 0.25, prior = "noninformative")

  expect_no_error(plot(model, dims = 1, n_grid = 8))
  expect_no_error(plot(model, dims = c(1, 2), n_grid = 8))
  expect_error(plot(model, only_mean = c(TRUE, FALSE)), "`only_mean` must be a single logical value")
  expect_error(plot(model, n_grid = 0), "'n_grid' must be a positive integer")
  expect_error(plot(model, dims = c(1.1, 2)), "`dims` must be an integer vector")
  expect_error(plot(model, dims = integer(0)), "`dims` must have length 1 or 2")
  expect_error(plot(model, dims = c(1, 1)), "`dims` cannot contain duplicate indices")
  expect_error(plot(model, dims = 3), "must be within the range")
})


test_that("plot.DKP supports ggplot engine for 1D/2D and validates engine", {
  set.seed(2026)
  skip_if_not_installed("ggplot2")

  X1 <- matrix(runif(60), ncol = 1)
  cl1 <- sample(1:3, nrow(X1), replace = TRUE)
  Y1 <- matrix(0, nrow(X1), 3)
  Y1[cbind(seq_len(nrow(X1)), cl1)] <- 1
  model1 <- fit_DKP(X1, Y1, prior = "noninformative", theta = 0.25)
  expect_no_error(plot(model1, n_grid = 8, engine = "ggplot"))

  X <- matrix(runif(80), ncol = 2)
  cl <- sample(1:3, nrow(X), replace = TRUE)
  Y <- matrix(0, nrow(X), 3)
  Y[cbind(seq_len(nrow(X)), cl)] <- 1
  model <- fit_DKP(X, Y, theta = 0.25, prior = "noninformative")

  expect_no_error(plot(model, n_grid = 8, engine = "ggplot"))
  expect_error(plot(model, engine = "bad_engine"),
               regexp = "one of.*base.*ggplot",
               ignore.case = TRUE)
})
