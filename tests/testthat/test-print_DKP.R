# This test file validates the functionality of the `print` S3 method for
# DKP objects and their related output classes (summary, predict, simulate).
# It captures console output so that print methods do not pollute test output.

test_that("print.DKP methods run without errors and produce expected output", {
  set.seed(123)

  # -------------------------------------------------------------------------
  # Setup: Create a DKP model and related objects
  # -------------------------------------------------------------------------

  true_pi_fun <- function(X) {
    p1 <- 1 / (1 + exp(-3 * X))
    p2 <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
    matrix(c(p1 / 2, p2 / 2, 1 - (p1 + p2) / 2), nrow = length(p1))
  }

  n <- 30
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- matrix(runif(n, Xbounds[1], Xbounds[2]), ncol = 1)
  true_pi <- true_pi_fun(X)
  m <- sample(150, n, replace = TRUE)

  Y <- t(sapply(seq_len(n), function(i) {
    rmultinom(1, size = m[i], prob = true_pi[i, ])
  }))

  model <- fit_DKP(
    X, Y,
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
    "Dirichlet Kernel Process"
  )

  expect_output(
    print(summary_model),
    "Dirichlet Kernel Process"
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


test_that("print.DKP handles new-data, multi-class and high-dimensional branches", {
  set.seed(100)

  X <- matrix(runif(45), ncol = 3)
  Y <- matrix(0, nrow(X), 4)

  cls <- sample(1:4, nrow(X), replace = TRUE)
  Y[cbind(seq_len(nrow(X)), cls)] <- 1

  model <- fit_DKP(
    X, Y,
    prior = "fixed",
    r0 = 6,
    p0 = rep(0.25, 4),
    theta = c(0.35, 0.4, 0.45),
    isotropic = FALSE
  )

  Xnew <- X[1:5, , drop = FALSE]

  expect_output(
    print(model),
    "Dirichlet Kernel Process"
  )

  expect_output(
    print(summary(model)),
    "Dirichlet Kernel Process"
  )

  expect_output(
    print(predict(model, Xnew = Xnew)),
    "Prediction results"
  )

  expect_output(
    print(simulate(model, Xnew = Xnew, nsim = 5)),
    "Simulation results"
  )
})
