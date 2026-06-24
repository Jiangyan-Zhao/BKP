# This test file validates the functionality of the `print` S3 method for
# BKP objects and their related output classes (summary, predict, simulate).
# It captures console output so that print methods do not pollute test output.

test_that("print.BKP methods run without errors and produce expected output", {
  set.seed(123)

  # -------------------------------------------------------------------------
  # Setup: Create a BKP model and related objects
  # -------------------------------------------------------------------------

  true_pi_fun <- function(x) {
    (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
  }

  n <- 30
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- matrix(runif(n, Xbounds[1], Xbounds[2]), ncol = 1)
  true_pi <- true_pi_fun(X)
  m <- sample(100, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = true_pi)

  model <- fit_BKP(
    X, y, m,
    Xbounds = Xbounds,
    prior = "noninformative"
  )

  summary_model <- summary(model)
  predict_model <- predict(model)
  simulate_model <- simulate(model)

  # -------------------------------------------------------------------------
  # Test Cases: Verify print methods while capturing console output
  # -------------------------------------------------------------------------

  expect_output(
    print(model),
    "Beta Kernel Process"
  )

  expect_output(
    print(summary_model),
    "Beta Kernel Process"
  )

  expect_output(
    print(predict_model),
    "Prediction results"
  )

  expect_output(
    print(simulate_model),
    "Simulation results"
  )
})


test_that("print.BKP handles new-data and high-dimensional branches", {
  set.seed(100)

  X <- matrix(runif(36), ncol = 3)
  m <- rep(1, nrow(X))
  y <- rbinom(nrow(X), size = 1, prob = 0.5)

  model <- fit_BKP(
    X, y, m,
    prior = "fixed",
    r0 = 5,
    p0 = 0.5,
    theta = c(0.3, 0.4, 0.5),
    isotropic = FALSE
  )

  Xnew <- X[1:5, , drop = FALSE]

  expect_output(
    print(model),
    "Beta Kernel Process"
  )

  expect_output(
    print(summary(model)),
    "Beta Kernel Process"
  )

  expect_output(
    print(predict(model, Xnew = Xnew, threshold = 0.4)),
    "Prediction results"
  )

  expect_output(
    print(simulate(model, Xnew = Xnew, nsim = 5, threshold = 0.4)),
    "Simulation results"
  )
})
