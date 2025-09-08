# This test file validates the functionality of the `print` S3 method for
# DKP objects and their related output classes (summary, predict, simulate).
# It ensures that these methods run without errors.

test_that("print.DKP methods run without errors", {
  # Set a seed for reproducibility
  set.seed(123)

  # -------------------------------------------------------------------------
  # Setup: Create a DKP model and related objects
  # -------------------------------------------------------------------------

  # Define true class probability function (3-class)
  true_pi_fun <- function(X) {
    p1 <- 1/(1+exp(-3*X))
    p2 <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
    return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
  }

  n <- 30
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- matrix(runif(n, Xbounds[1], Xbounds[2]), ncol = 1)
  true_pi <- true_pi_fun(X)
  m <- sample(150, n, replace = TRUE)
  Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))

  # Fit DKP model
  model <- fit_DKP(X, Y, Xbounds = Xbounds, prior = "noninformative")

  # Generate summary, predict, and simulate objects
  summary_model <- summary(model)
  predict_model <- predict(model)
  simulate_model <- simulate(model)

  # -------------------------------------------------------------------------
  # Test Cases: Verify print methods
  # -------------------------------------------------------------------------

  # Test print.DKP
  expect_no_error(print(model))

  # Test print.summary_DKP
  expect_no_error(print(summary_model))

  # Test print.predict_DKP
  expect_no_error(print(predict_model))

  # Test print.simulate_DKP
  expect_no_error(print(simulate_model))
})
