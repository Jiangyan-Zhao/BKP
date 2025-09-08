# This test file validates the functionality of the `quantile` S3 method for
# BKP objects. It checks for correct output and proper input validation.

test_that("quantile.BKP returns correct posterior quantiles", {
  # Set a seed for reproducibility
  set.seed(123)

  # -------------------------------------------------------------------------
  # Setup: Create a real BKP model object using fit_BKP
  # -------------------------------------------------------------------------

  # Define true success probability function
  true_pi_fun <- function(x) {
    (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
  }

  n <- 30
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- tgp::lhs(n = n, rect = Xbounds)
  true_pi <- true_pi_fun(X)
  m <- sample(100, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = true_pi)

  # Fit BKP model (this will be the object to test)
  model <- fit_BKP(X, y, m, Xbounds = Xbounds)

  # -------------------------------------------------------------------------
  # Test Cases: Verify quantile function
  # -------------------------------------------------------------------------

  # Get posterior parameters from the fitted model
  alpha_n <- model$alpha_n
  beta_n <- model$beta_n

  # Test with default quantiles (0.025, 0.5, 0.975)
  quantiles_default <- quantile(model)
  # Use sapply to ensure the expected output has the correct matrix dimensions
  expected_default <- sapply(c(0.025, 0.5, 0.975), function(p) qbeta(p, alpha_n, beta_n))
  dimnames(expected_default) <- dimnames(quantiles_default)

  expect_equal(quantiles_default, expected_default)
  expect_true(is.matrix(quantiles_default))
  expect_equal(dim(quantiles_default), c(n, 3))

  # Test with a single quantile
  quantiles_single <- quantile(model, probs = 0.5)
  expected_single <- qbeta(0.5, alpha_n, beta_n)

  expect_equal(quantiles_single, as.vector(expected_single))
  expect_true(is.numeric(quantiles_single))
  expect_equal(length(quantiles_single), n)

  # Test with a custom vector of quantiles
  custom_probs <- c(0.1, 0.9)
  quantiles_custom <- quantile(model, probs = custom_probs)
  # Use sapply to ensure the expected output has the correct matrix dimensions
  expected_custom <- sapply(custom_probs, function(p) qbeta(p, alpha_n, beta_n))
  dimnames(expected_custom) <- dimnames(quantiles_custom)

  expect_equal(quantiles_custom, expected_custom)
  expect_true(is.matrix(quantiles_custom))
  expect_equal(dim(quantiles_custom), c(n, 2))
})

test_that("quantile.BKP handles input validation correctly", {
  # Create a dummy BKP model object
  alpha_n <- rep(1, 5)
  beta_n <- rep(1, 5)
  mock_model <- list(alpha_n = alpha_n, beta_n = beta_n)
  class(mock_model) <- "BKP"

  # Test for invalid probs values
  expect_error(quantile(mock_model, probs = c(0.5, 1.1)),
               "'probs' must be a numeric vector with all values in \\[0, 1\\].")
  expect_error(quantile(mock_model, probs = -0.1),
               "'probs' must be a numeric vector with all values in \\[0, 1\\].")
  expect_error(quantile(mock_model, probs = "a"),
               "'probs' must be a numeric vector with all values in \\[0, 1\\].")
})
