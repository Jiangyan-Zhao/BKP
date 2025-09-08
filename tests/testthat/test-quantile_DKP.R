# This test file validates the functionality of the `quantile` S3 method for
# DKP objects. It checks for correct output and proper input validation,
# using a real DKP model object.

test_that("quantile.DKP returns correct posterior quantiles", {
  # Set a seed for reproducibility
  set.seed(123)

  # -------------------------------------------------------------------------
  # Setup: Create a real DKP model object using fit_DKP
  # -------------------------------------------------------------------------

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
  Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))

  # Fit DKP model (this will be the object to test)
  model <- fit_DKP(X, Y, Xbounds = Xbounds)

  # -------------------------------------------------------------------------
  # Test Cases: Verify quantile function
  # -------------------------------------------------------------------------

  # Get posterior parameters from the fitted model
  alpha_n <- model$alpha_n
  row_sum <- rowSums(alpha_n)

  # Test with default quantiles (0.025, 0.5, 0.975)
  quantiles_default <- quantile(model)
  # Use sapply to ensure the expected output has the correct dimensions
  expected_default <- array(NA, dim = dim(quantiles_default),
                            dimnames = dimnames(quantiles_default))
  for(j in 1:ncol(alpha_n)) {
    expected_default[, j, ] <- t(sapply(1:n, function(i) {
      qbeta(c(0.025, 0.5, 0.975), alpha_n[i, j], row_sum[i] - alpha_n[i, j])
    }))
  }

  expect_equal(quantiles_default, expected_default)
  expect_true(is.array(quantiles_default))
  expect_equal(dim(quantiles_default), c(n, 3, 3))

  # Test with a single quantile
  quantiles_single <- quantile(model, probs = 0.5)
  expected_single <- matrix(NA, nrow = n, ncol = ncol(alpha_n))
  for(j in 1:ncol(alpha_n)) {
    expected_single[, j] <- sapply(1:n, function(i) {
      qbeta(0.5, alpha_n[i, j], row_sum[i] - alpha_n[i, j])
    })
  }

  # Manually clear the dimnames to allow for proper comparison
  dimnames(quantiles_single) <- NULL
  dimnames(expected_single) <- NULL

  expect_equal(quantiles_single, expected_single)
  expect_true(is.matrix(quantiles_single))
  expect_equal(dim(quantiles_single), c(n, 3))

  # Test with a custom vector of quantiles
  custom_probs <- c(0.1, 0.9)
  quantiles_custom <- quantile(model, probs = custom_probs)
  expected_custom <- array(NA, dim = dim(quantiles_custom),
                           dimnames = dimnames(quantiles_custom))
  for(j in 1:ncol(alpha_n)) {
    expected_custom[, j, ] <- t(sapply(1:n, function(i) {
      qbeta(custom_probs, alpha_n[i, j], row_sum[i] - alpha_n[i, j])
    }))
  }

  expect_equal(quantiles_custom, expected_custom)
  expect_true(is.array(quantiles_custom))
  expect_equal(dim(quantiles_custom), c(n, 3, 2))
})

test_that("quantile.DKP handles input validation correctly", {
  # Set a seed for reproducibility
  set.seed(123)

  # Create a real DKP model object
  # (This ensures the test is performed on a genuine object)
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
  Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))

  model <- fit_DKP(X, Y, Xbounds = Xbounds)

  # Test for invalid probs values on the real model. The test should PASS
  # because the function correctly throws an error.
  # We use a regex for a more robust check on the error message.
  expect_error(quantile(model, probs = c(0.5, 1.1)),
               "'probs' must be a numeric vector with all values in \\[0, 1\\].")
  expect_error(quantile(model, probs = -0.1),
               "'probs' must be a numeric vector with all values in \\[0, 1\\].")
  expect_error(quantile(model, probs = "a"),
               "'probs' must be a numeric vector with all values in \\[0, 1\\].")
})
