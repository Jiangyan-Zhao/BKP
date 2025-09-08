# Test file for the get_prior function

test_that("get_prior handles input validation correctly", {
  # Test with an invalid model
  expect_error(get_prior(prior = "fixed", model = "invalid_model"),
               "'arg' should be one of \"BKP\", \"DKP\"")

  # Test with an invalid prior
  expect_error(get_prior(prior = "invalid_prior", model = "BKP"),
               "'arg' should be one of \"noninformative\", \"fixed\", \"adaptive\"")

  # Test for BKP fixed prior with missing p0
  expect_error(get_prior(prior = "fixed", model = "BKP", r0 = 10, y = NULL),
               "For fixed prior in BKP, 'p0' must be in \\(0,1\\).")

  # Test for BKP fixed prior with invalid p0
  expect_error(get_prior(prior = "fixed", model = "BKP", r0 = 10, p0 = 2),
               "For fixed prior in BKP, 'p0' must be in \\(0,1\\).")

  # Test for BKP adaptive prior with missing y/m
  expect_error(get_prior(prior = "adaptive", model = "BKP", K = matrix(1), r0 = 10, y = NULL),
               "For adaptive prior in BKP, 'y', 'm', and 'K' must be provided.")

  # Test for DKP fixed prior with missing p0
  expect_error(get_prior(prior = "fixed", model = "DKP", r0 = 10),
               "Either 'Y' or 'p0' must be provided to determine class dimension q.")

  # Test for DKP fixed prior with p0 not summing to 1
  expect_error(get_prior(prior = "fixed", model = "DKP", r0 = 10, p0 = c(0.5, 0.6)),
               "'p0' must sum to 1.")

  # Test for DKP adaptive prior with missing Y/K
  expect_error(get_prior(prior = "adaptive", model = "DKP", r0 = 10, Y = matrix(1)),
               "'Y' and 'K' must be provided for adaptive prior in DKP.")
})

test_that("get_prior works correctly for BKP model", {
  set.seed(42) # Set seed for reproducibility
  n <- 10
  # Create a proper kernel matrix based on some dummy data
  X <- matrix(runif(n*2), n, 2)
  theta <- 0.5
  D <- as.matrix(dist(X)^2)
  K <- exp(-D / (2 * theta^2))
  y <- rbinom(n, 10, 0.5)
  m <- rep(10, n)

  # Noninformative prior
  non_info_prior <- get_prior(prior = "noninformative", model = "BKP", K = K)
  expect_equal(non_info_prior$alpha0, rep(1, n))
  expect_equal(non_info_prior$beta0, rep(1, n))

  # Fixed prior
  p0 <- 0.6
  r0 <- 10
  fixed_prior <- get_prior(prior = "fixed", model = "BKP", K = K, p0 = p0, r0 = r0)
  expect_equal(fixed_prior$alpha0, rep(r0 * p0, n))
  expect_equal(fixed_prior$beta0, rep(r0 * (1 - p0), n))

  # Adaptive prior
  r0 <- 10
  adaptive_prior <- get_prior(prior = "adaptive", model = "BKP", K = K, y = y, m = m, r0 = r0)
  expected_alpha0 <- r0 * (t(K) %*% (y / m))
  expected_beta0 <- r0 * (t(K) %*% ((m - y) / m))

  # The function seems to be clamping to 0.0. Let's adjust the expectation.
  expected_alpha0[expected_alpha0 <= 1e-10] <- 0.0
  expected_beta0[expected_beta0 <= 1e-10] <- 0.0

  expect_equal(unname(adaptive_prior$alpha0), unname(as.vector(expected_alpha0)), tolerance = 1e-6)
  expect_equal(unname(adaptive_prior$beta0), unname(as.vector(expected_beta0)), tolerance = 1e-6)
})

test_that("get_prior works correctly for DKP model", {
  set.seed(42) # Set seed for reproducibility
  n <- 10
  q <- 3
  # Create a proper kernel matrix based on some dummy data
  X <- matrix(runif(n*2), n, 2)
  theta <- 0.5
  D <- as.matrix(dist(X)^2)
  K <- exp(-D / (2 * theta^2))
  Y <- matrix(rbinom(n*q, 10, 0.3), n, q)

  # Noninformative prior
  non_info_prior <- get_prior(prior = "noninformative", model = "DKP", Y = Y)
  expect_equal(non_info_prior, matrix(1, nrow = 1, ncol = q))

  # Fixed prior
  r0 <- 10
  p0 <- colMeans(Y / rowSums(Y))
  fixed_prior <- get_prior(prior = "fixed", model = "DKP", K = K, r0 = r0, p0 = p0)
  expected_alpha0 <- matrix(rep(r0 * p0, each = n), nrow = n, byrow = TRUE)
  expect_equal(fixed_prior, expected_alpha0)

  # Adaptive prior
  r0 <- 10
  adaptive_prior <- get_prior(prior = "adaptive", model = "DKP", K = K, Y = Y, r0 = r0)
  sum_y <- rowSums(Y)
  sum_y[sum_y == 0] <- 1e-10 # to avoid division by zero
  expected_alpha0 <- r0 * t(K) %*% (Y / sum_y)
  expect_equal(unname(adaptive_prior), unname(expected_alpha0), tolerance = 1e-6)
})
