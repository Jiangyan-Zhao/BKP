# This test file validates the functionality of the `simulate` S3 method for
# BKP objects. It checks for correct output dimensions, value ranges, and
# proper input validation.

test_that("simulate.BKP returns correct posterior simulations and dimensions", {
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
  # Test Cases
  # -------------------------------------------------------------------------

  # Case 1: Basic simulation on training data
  nsim_test <- 50
  sim_result <- simulate(model, nsim = nsim_test)

  # Check dimensions of samples
  expect_equal(dim(sim_result$samples), c(n, nsim_test))

  # Check if samples are within [0, 1]
  expect_true(all(sim_result$samples >= 0 & sim_result$samples <= 1))

  # Check for `class` output when threshold is NULL
  expect_null(sim_result$class)


  # Case 2: Simulation on new data
  n_new <- 15
  Xnew <- tgp::lhs(n = n_new, rect = Xbounds)
  sim_new_result <- simulate(model, Xnew = Xnew, nsim = 10)

  # Check dimensions of samples for new data
  expect_equal(dim(sim_new_result$samples), c(n_new, 10))

  # Check for `class` output when threshold is NULL
  expect_null(sim_new_result$class)


  # Case 3: Simulation with classification threshold
  sim_class_result <- simulate(model, nsim = 20, threshold = 0.5)

  # Check if `class` output exists and has correct dimensions
  expect_true(!is.null(sim_class_result$class))
  expect_equal(dim(sim_class_result$class), c(n, 20))

  # Check if classification values are binary (0 or 1)
  expect_true(all(sim_class_result$class %in% c(0, 1)))

  # Check if classification logic holds
  samples_matrix <- sim_class_result$samples
  class_matrix <- sim_class_result$class
  expect_true(all(class_matrix == ifelse(samples_matrix > 0.5, 1, 0)))


  # Case 4: Test for reproducibility with a fixed seed
  sim_seed1 <- simulate(model, nsim = 10, seed = 42)
  sim_seed2 <- simulate(model, nsim = 10, seed = 42)
  expect_equal(sim_seed1$samples, sim_seed2$samples)

  # Test that different seeds produce different results
  sim_seed3 <- simulate(model, nsim = 10, seed = 1)
  expect_false(isTRUE(all.equal(sim_seed1$samples, sim_seed3$samples)))
})

test_that("simulate.BKP handles input validation correctly", {
  # Set up a test model
  set.seed(123)
  true_pi_fun <- function(x) (1 + exp(-x^2)) / 2
  n <- 10
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- tgp::lhs(n = n, rect = Xbounds)
  true_pi <- true_pi_fun(X)
  m <- sample(100, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = true_pi)
  model <- fit_BKP(X, y, m, Xbounds = Xbounds)

  # Case 5: Input validation tests
  # nsim must be a positive integer
  expect_error(simulate(model, nsim = 0), "`nsim` must be a positive integer.")
  expect_error(simulate(model, nsim = -10), "`nsim` must be a positive integer.")
  expect_error(simulate(model, nsim = 5.5), "`nsim` must be a positive integer.")

  # threshold must be a single numeric value in (0, 1)
  expect_error(simulate(model, threshold = 1.1), "`threshold` must be a numeric value strictly between 0 and 1 (e.g., 0.5).", fixed = TRUE)
  expect_error(simulate(model, threshold = c(0.1, 0.9)), "`threshold` must be a numeric value strictly between 0 and 1 (e.g., 0.5).", fixed = TRUE)

  # Xnew must be a matrix
  expect_error(simulate(model, Xnew = "invalid"), "'Xnew' must be numeric.")
})
