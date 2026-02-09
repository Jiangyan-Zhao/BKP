# This test file validates the functionality of the `print` S3 method for
# BKP objects and their related output classes (summary, predict, simulate).
# It ensures that these methods run without errors.

test_that("print.BKP methods run without errors", {
  # Set a seed for reproducibility
  set.seed(123)

  # -------------------------------------------------------------------------
  # Setup: Create a BKP model and related objects
  # -------------------------------------------------------------------------

  # Define true success probability function
  true_pi_fun <- function(x) {
    (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
  }

  n <- 30
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- matrix(runif(n, Xbounds[1], Xbounds[2]), ncol = 1)
  true_pi <- true_pi_fun(X)
  m <- sample(100, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = true_pi)

  # Fit BKP model
  model <- fit_BKP(X, y, m, Xbounds = Xbounds, prior = "noninformative")

  # Generate summary, predict, and simulate objects
  summary_model <- summary(model)
  predict_model <- predict(model)
  simulate_model <- simulate(model)

  # -------------------------------------------------------------------------
  # Test Cases: Verify print methods
  # -------------------------------------------------------------------------

  # Test print.BKP
  expect_no_error(print(model))

  # Test print.summary_BKP
  expect_no_error(print(summary_model))

  # Test print.predict_BKP
  expect_no_error(print(predict_model))

  # Test print.simulate_BKP
  expect_no_error(print(simulate_model))
})

test_that("print.BKP handles new-data and high-dimensional branches", {
  set.seed(100)

  X <- matrix(runif(36), ncol = 3)
  m <- rep(1, nrow(X))
  y <- rbinom(nrow(X), size = 1, prob = 0.5)
  model <- fit_BKP(X, y, m, prior = "fixed", r0 = 5, p0 = 0.5, theta = c(0.3, 0.4, 0.5), isotropic = FALSE)

  expect_no_error(print(model))
  expect_no_error(print(summary(model)))
  expect_no_error(print(predict(model, Xnew = X[1:5, , drop = FALSE], threshold = 0.4)))
  expect_no_error(print(simulate(model, Xnew = X[1:5, , drop = FALSE], nsim = 5, threshold = 0.4)))
})
