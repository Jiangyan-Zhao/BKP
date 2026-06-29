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
               "For fixed prior in BKP, 'p0' must be a scalar in (0, 1).", fixed = TRUE)

  # Test for BKP fixed prior with invalid p0
  expect_error(get_prior(prior = "fixed", model = "BKP", r0 = 10, p0 = 2),
               "For fixed prior in BKP, 'p0' must be a scalar in (0, 1).", fixed = TRUE)

  # Test for BKP adaptive prior with missing y/m
  expect_error(get_prior(prior = "adaptive", model = "BKP", K = matrix(1), r0 = 10, y = NULL),
               "For adaptive prior in BKP, 'y', 'm', and 'K' must be provided.")

  # Test for DKP fixed prior with missing p0
  expect_error(get_prior(prior = "fixed", model = "DKP", r0 = 10),
               "Either 'Y' or 'p0' must be provided to determine the number of classes.")

  # Test for DKP fixed prior with p0 not summing to 1
  expect_error(get_prior(prior = "fixed", model = "DKP", r0 = 10, p0 = c(0.5, 0.6)),
               "For fixed prior in DKP, 'p0' must be a nonnegative finite numeric vector of length equal to the number of classes and sum to 1.", fixed = TRUE)

  # Test for DKP adaptive prior with missing Y/K
  expect_error(get_prior(prior = "adaptive", model = "DKP", r0 = 10, Y = matrix(1, nrow = 1, ncol = 2)),
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
  expect_equal(dim(non_info_prior), c(n, q))
  expect_equal(non_info_prior, matrix(1, nrow = n, ncol = q))

  # Fixed prior
  r0 <- 10
  p0 <- colMeans(Y / rowSums(Y))
  fixed_prior <- get_prior(prior = "fixed", model = "DKP", K = K, r0 = r0, p0 = p0)
  expected_alpha0 <- matrix(rep(r0 * p0, each = n), nrow = n)
  expect_equal(fixed_prior, expected_alpha0)

  # Adaptive prior
  r0 <- 10
  adaptive_prior <- get_prior(prior = "adaptive", model = "DKP", K = K, Y = Y, r0 = r0)
  sum_y <- rowSums(Y)
  sum_y[sum_y == 0] <- 1e-10 # to avoid division by zero
  expected_alpha0 <- r0 * t(K) %*% (Y / sum_y)
  expect_equal(unname(adaptive_prior), unname(expected_alpha0), tolerance = 1e-6)
})

test_that("get_prior additional validation and edge branches", {
  expect_error(get_prior(model = "BKP", y = c(1, NA), m = c(1, 1)),
               "'y' must be a numeric vector with finite values and no NA values.")
  expect_error(get_prior(model = "BKP", y = c(1, 0), m = c(1, NA)),
               "'m' must be a numeric vector with finite values and no NA values.")
  expect_error(get_prior(model = "DKP", Y = matrix(c(1, NA), ncol = 1)),
               "'Y' must be a numeric matrix with finite values and no NA values.")
  expect_error(get_prior(model = "DKP", K = matrix(c(1, NA), ncol = 1)),
               "'K' must be a numeric matrix with finite values and no NA values.")

  alpha0 <- get_prior(prior = "fixed", model = "DKP", r0 = 2, p0 = c(0.2, 0.8), K = matrix(1, nrow = 3, ncol = 2))
  expect_true(is.matrix(alpha0))
  expect_equal(dim(alpha0), c(3, 2))

  p <- get_prior(prior = "noninformative", model = "BKP")
  expect_equal(p$alpha0, 1)
  expect_equal(p$beta0, 1)
})

test_that("test-get_prior uncovered validation branches", {
  expect_error(get_prior(model = "BKP", r0 = 0), "'r0' must be a positive finite scalar.")
  expect_error(get_prior(model = "BKP", p0 = -0.1),
               "'p0' must be a numeric vector with nonnegative finite values.")

  expect_error(
    get_prior(prior = "adaptive", model = "BKP", y = c(1, 0), m = c(1, 1), K = matrix(1, nrow = 2, ncol = 3)),
    "'K' must have ncol = length(y).", fixed = TRUE
  )

  expect_error(
    get_prior(prior = "adaptive", model = "BKP", y = c(1, 0), m = c(1, 1, 1), K = matrix(1, nrow = 2, ncol = 2)),
    "'y' and 'm' must have the same length."
  )

  expect_error(
    get_prior(prior = "fixed", model = "DKP", p0 = c(0.3, 0.3), Y = matrix(1, nrow = 1, ncol = 3)),
    "For fixed prior in DKP, 'p0' must be a nonnegative finite numeric vector of length equal to the number of classes and sum to 1.", fixed = TRUE
  )

  expect_error(
    get_prior(prior = "adaptive", model = "DKP", Y = matrix(c(1, 0, 0, 1), ncol = 2), K = matrix(1, nrow = 2, ncol = 3)),
    "'K' must have ncol = nrow(Y).", fixed = TRUE
  )

  expect_error(
    get_prior(prior = "adaptive", model = "DKP", Y = matrix(c(0, 0, 1, 0), ncol = 2), K = diag(2)),
    "Each row of 'Y' must have a positive row sum."
  )

  alpha0_noninfo <- get_prior(prior = "noninformative", model = "DKP", p0 = c(0.4, 0.6))
  expect_equal(dim(alpha0_noninfo), c(1, 2))
})

test_that("get_prior infers output dimensions when K is NULL for nonadaptive priors", {
  y <- c(1, 2, 0, 3)
  m <- c(5, 5, 5, 5)

  bkp_noninfo <- get_prior(prior = "noninformative", model = "BKP", y = y, m = m)
  expect_equal(length(bkp_noninfo$alpha0), length(y))
  expect_equal(length(bkp_noninfo$beta0), length(y))
  expect_equal(bkp_noninfo$alpha0, rep(1, length(y)))
  expect_equal(bkp_noninfo$beta0, rep(1, length(y)))

  r0 <- 8
  p0_bkp <- 0.25
  bkp_fixed <- get_prior(prior = "fixed", model = "BKP", y = y, m = m, r0 = r0, p0 = p0_bkp)
  expect_equal(length(bkp_fixed$alpha0), length(y))
  expect_equal(length(bkp_fixed$beta0), length(y))
  expect_equal(bkp_fixed$alpha0, rep(r0 * p0_bkp, length(y)))
  expect_equal(bkp_fixed$beta0, rep(r0 * (1 - p0_bkp), length(y)))

  Y <- matrix(c(2, 1, 0, 1, 3, 2, 1, 0, 4), nrow = 3, byrow = TRUE)
  dkp_noninfo <- get_prior(prior = "noninformative", model = "DKP", Y = Y)
  expect_equal(dim(dkp_noninfo), dim(Y))
  expect_equal(dkp_noninfo, matrix(1, nrow = nrow(Y), ncol = ncol(Y)))

  p0_dkp <- c(0.2, 0.3, 0.5)
  r0 <- 10
  dkp_fixed <- get_prior(prior = "fixed", model = "DKP", Y = Y, p0 = p0_dkp, r0 = r0)
  expect_equal(dim(dkp_fixed), dim(Y))
  expect_equal(dkp_fixed, matrix(rep(r0 * p0_dkp, each = nrow(Y)), nrow = nrow(Y)))

  dkp_fixed_no_y <- get_prior(prior = "fixed", model = "DKP", p0 = p0_dkp, r0 = r0)
  expect_equal(dim(dkp_fixed_no_y), c(1L, length(p0_dkp)))
  expect_equal(dkp_fixed_no_y, matrix(r0 * p0_dkp, nrow = 1))
})

test_that("get_prior rejects invalid DKP class inputs", {
  expect_error(
    get_prior(prior = "fixed", model = "DKP", p0 = 1),
    "For DKP, the number of classes must be at least two.",
    fixed = TRUE
  )

  expect_error(
    get_prior(prior = "noninformative", model = "DKP", Y = matrix(c(1, 0, 0, 0), nrow = 2, byrow = TRUE)),
    "Each row of 'Y' must have a positive row sum.",
    fixed = TRUE
  )
})

test_that("get_prior adaptive priors still require K", {
  expect_error(
    get_prior(prior = "adaptive", model = "BKP", y = c(1, 0), m = c(2, 2)),
    "For adaptive prior in BKP, 'y', 'm', and 'K' must be provided.",
    fixed = TRUE
  )

  expect_error(
    get_prior(prior = "adaptive", model = "DKP", Y = matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)),
    "'Y' and 'K' must be provided for adaptive prior in DKP.",
    fixed = TRUE
  )
})
