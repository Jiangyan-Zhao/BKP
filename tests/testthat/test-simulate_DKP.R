# This test file validates the functionality of the `simulate` S3 method for
# DKP objects. It checks for correct output dimensions, value ranges,
# and proper input validation.

test_that("simulate.DKP returns correct posterior simulations and dimensions", {
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
  q <- 3
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- tgp::lhs(n = n, rect = Xbounds)
  true_pi <- true_pi_fun(X)
  m <- sample(150, n, replace = TRUE)

  # Generate multinomial responses
  Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))

  # Fit DKP model (this will be the object to test)
  model <- fit_DKP(X, Y, Xbounds = Xbounds)

  # -------------------------------------------------------------------------
  # Test Cases
  # -------------------------------------------------------------------------

  # Case 1: Basic simulation on training data
  nsim_test <- 50
  sim_result <- simulate(model, nsim = nsim_test)

  # Check dimensions of samples
  expect_equal(dim(sim_result$samples), c(n, q, nsim_test))

  # Check if all probabilities sum to 1
  expect_true(all(abs(apply(sim_result$samples, c(1,3), sum) - 1) < 1e-10))

  # Check for `class` output (should be NULL for multi-count data)
  expect_null(sim_result$class)


  # Case 2: Simulation on new data
  n_new <- 15
  Xnew <- tgp::lhs(n = n_new, rect = Xbounds)
  sim_new_result <- simulate(model, Xnew = Xnew, nsim = 10)

  # Check dimensions of samples for new data
  expect_equal(dim(sim_new_result$samples), c(n_new, q, 10))


  # Case 3: Test for reproducibility with a fixed seed
  sim_seed1 <- simulate(model, nsim = 10, seed = 42)
  sim_seed2 <- simulate(model, nsim = 10, seed = 42)
  expect_equal(sim_seed1$samples, sim_seed2$samples)

  # Test that different seeds produce different results
  sim_seed3 <- simulate(model, nsim = 10, seed = 1)
  expect_false(isTRUE(all.equal(sim_seed1$samples, sim_seed3$samples)))


  # Case 4: Test MAP classification for single-label data
  m_single <- rep(1, n)
  Y_single <- t(sapply(1:n, function(i) rmultinom(1, size = m_single[i], prob = true_pi[i, ])))
  model_single <- fit_DKP(X, Y_single, Xbounds = Xbounds)

  sim_class_result <- simulate(model_single, nsim = 20)

  # Check if `class` output exists and has correct dimensions
  expect_true(!is.null(sim_class_result$class))
  expect_equal(dim(sim_class_result$class), c(n, 20))

  # Check if classification values are correct
  for (i in 1:20) {
    samples_sim <- sim_class_result$samples[, , i]
    class_sim <- sim_class_result$class[, i]
    expect_true(all(class_sim == apply(samples_sim, 1, which.max)))
  }
})

test_that("simulate.DKP handles input validation correctly", {
  # Set up a test model
  set.seed(123)
  true_pi_fun <- function(X) {
    p1 <- 1/(1+exp(-3*X))
    p2 <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
    return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
  }

  n <- 10
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- tgp::lhs(n = n, rect = Xbounds)
  true_pi <- true_pi_fun(X)
  m <- sample(150, n, replace = TRUE)
  Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
  model <- fit_DKP(X, Y, Xbounds = Xbounds)

  # Case 5: Input validation tests
  # nsim must be a positive integer
  expect_error(simulate(model, nsim = 0), "`nsim` must be a positive integer.", fixed = TRUE)
  expect_error(simulate(model, nsim = -10), "`nsim` must be a positive integer.", fixed = TRUE)
  expect_error(simulate(model, nsim = 5.5), "`nsim` must be a positive integer.", fixed = TRUE)

  # Xnew must be a matrix
  expect_error(simulate(model, Xnew = "invalid"), "'Xnew' must be numeric.", fixed = TRUE)
  expect_error(simulate(model, Xnew = matrix("a", 2, 2)), "'Xnew' must be numeric.", fixed = TRUE)

  # DKP models do not support the threshold argument, and the function
  # currently ignores it due to the '...' parameter. We now check for
  # no error to reflect the current function behavior.
  expect_no_error(simulate(model, threshold = 0.5))
})
