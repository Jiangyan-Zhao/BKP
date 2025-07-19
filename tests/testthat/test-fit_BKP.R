# -------------------------- Test Setup ---------------------------
set.seed(123)
n_test <- 50
d_test <- 2
X_test <- matrix(runif(n_test * d_test, -5, 5), n_test, d_test)
y_test <- sample(0:10, n_test, replace = TRUE)
m_test <- rep(10, n_test)

# Define a simple function that simulates a true probability surface for testing
true_prob_fun <- function(X) {
  # Example: 2D sine wave probability
  sin(X[,1]) * cos(X[,2]) / 2 + 0.5 # Scale to [0,1]
}

# -------------------------- Test Context: Basic Functionality ---------------------------

test_that("Basic Functionality: fit.BKP returns a 'BKP' class object with expected elements", {
  model <- fit.BKP(X_test, y_test, m_test)

  expect_s3_class(model, "BKP")
  expect_type(model, "list")

  expect_named(model, c("theta_opt", "kernel", "loss", "loss_min",
                        "X", "Xnorm", "Xbounds", "y", "m",
                        "prior", "r0", "p0", "alpha0", "beta0",
                        "alpha_n", "beta_n"))

  expect_true(all(model$X == X_test))
  expect_true(all(model$y == y_test))
  expect_true(all(model$m == m_test))
  expect_equal(length(model$theta_opt), d_test)
  expect_type(model$loss_min, "double")
})

test_that("Basic Functionality: fit.BKP works with different 'prior' types", {
  model_noninformative <- fit.BKP(X_test, y_test, m_test, prior = "noninformative")
  expect_equal(model_noninformative$prior, "noninformative")

  model_fixed <- fit.BKP(X_test, y_test, m_test, prior = "fixed", r0 = 5, p0 = 0.7)
  expect_equal(model_fixed$prior, "fixed")
  expect_equal(model_fixed$r0, 5)
  expect_equal(model_fixed$p0, 0.7)

  model_adaptive <- fit.BKP(X_test, y_test, m_test, prior = "adaptive", r0 = 10)
  expect_equal(model_adaptive$prior, "adaptive")
  expect_equal(model_adaptive$r0, 10)
})

test_that("Basic Functionality: fit.BKP works with different 'kernel' types", {
  model_gaussian <- fit.BKP(X_test, y_test, m_test, kernel = "gaussian")
  expect_equal(model_gaussian$kernel, "gaussian")

  model_matern52 <- fit.BKP(X_test, y_test, m_test, kernel = "matern52")
  expect_equal(model_matern52$kernel, "matern52")
})

test_that("Basic Functionality: fit.BKP works with different 'loss' functions", {
  model_brier <- fit.BKP(X_test, y_test, m_test, loss = "brier")
  expect_equal(model_brier$loss, "brier")

  model_logloss <- fit.BKP(X_test, y_test, m_test, loss = "log_loss")
  expect_equal(model_logloss$loss, "log_loss")
})

# -------------------------- Test Context: Input Validation ---------------------------

test_that("Input Validation: 'Xbounds' correctly normalizes inputs", {
  Xbounds_custom <- matrix(c(-10, 0, 10, 20), nrow = 2, byrow = TRUE)
  X_custom <- matrix(c(-5, 5, 0, 15), nrow = 2, byrow = TRUE) # Example values within custom bounds
  y_custom <- c(1,1)
  m_custom <- c(2,2)
  model <- fit.BKP(X_custom, y_custom, m_custom, Xbounds = Xbounds_custom)

  # Manually calculate expected normalized values
  Xnorm_expected <- (X_custom - matrix(Xbounds_custom[,1], nrow = nrow(X_custom), ncol = ncol(X_custom), byrow = TRUE)) /
    matrix((Xbounds_custom[,2] - Xbounds_custom[,1]), nrow = nrow(X_custom), ncol = ncol(X_custom), byrow = TRUE)

  expect_equal(model$Xnorm, Xnorm_expected)
})

# -------------------------- Test Context: Default Values ---------------------------
test_that("Default Values: Default 'prior' is 'noninformative'", {
  model <- fit.BKP(X_test, y_test, m_test)
  expect_equal(model$prior, "noninformative")
})

test_that("Default Values: Default 'r0' and 'p0' are used when prior is 'fixed'", {
  # r0=2, p0=0.5 are defaults in function signature
  model <- fit.BKP(X_test, y_test, m_test, prior = "fixed")
  expect_equal(model$r0, 2)
  expect_equal(model$p0, 0.5)
})

test_that("Default Values: Default 'kernel' is 'gaussian'", {
  model <- fit.BKP(X_test, y_test, m_test)
  expect_equal(model$kernel, "gaussian")
})

test_that("Default Values: Default 'loss' is 'brier'", {
  model <- fit.BKP(X_test, y_test, m_test)
  expect_equal(model$loss, "brier")
})



# -------------------------- Test Context: Edge Cases ---------------------------

test_that("Edge Cases: fit.BKP handles 1D input (d=1)", {
  X_1d <- matrix(runif(n_test, 0, 1), ncol = 1)
  model <- fit.BKP(X_1d, y_test, m_test)
  expect_s3_class(model, "BKP")
  expect_equal(ncol(model$X), 1)
  expect_equal(length(model$theta_opt), 1)
})

test_that("Edge Cases: fit.BKP handles small number of observations (n)", {
  X_small <- matrix(runif(5 * d_test), 5, d_test)
  y_small <- sample(0:10, 5, replace = TRUE)
  m_small <- rep(10, 5)
  model <- fit.BKP(X_small, y_small, m_small)
  expect_s3_class(model, "BKP")
  expect_equal(nrow(model$X), 5)
})

