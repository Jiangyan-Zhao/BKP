# -------------------------- Test Setup ---------------------------
set.seed(123)
# Define true success probability function
true_pi_fun <- function(x) {
  (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
}
d_test <- 1
n <- 30
Xbounds <- matrix(c(-2,2), nrow=1)
X <- tgp::lhs(n = n, rect = Xbounds)
true_pi <- true_pi_fun(X)
m <- sample(100, n, replace = TRUE)
y <- rbinom(n, size = m, prob = true_pi)

# -------------------------- Test Context: Basic Functionality ---------------------------

test_that("Basic Functionality: fit.BKP returns a 'BKP' class object with expected elements", {
  model <- fit.BKP(X, y, m, Xbounds)

  expect_s3_class(model, "BKP")
  expect_type(model, "list")

  expect_named(model, c("theta_opt", "kernel", "loss", "loss_min",
                        "X", "Xnorm", "Xbounds", "y", "m",
                        "prior", "r0", "p0", "alpha0", "beta0",
                        "alpha_n", "beta_n"))

  expect_equal(length(model$theta_opt), d_test)
  expect_type(model$loss_min, "double")
})

test_that("Basic Functionality: fit.BKP works with different 'prior' types", {
  model_noninformative <- fit.BKP(X, y, m, Xbounds, prior = "noninformative")
  expect_equal(model_noninformative$prior, "noninformative")

  model_fixed <- fit.BKP(X, y, m, Xbounds, prior = "fixed", r0 = 5, p0 = 0.7)
  expect_equal(model_fixed$prior, "fixed")
  expect_equal(model_fixed$r0, 5)
  expect_equal(model_fixed$p0, 0.7)

  model_adaptive <- fit.BKP(X, y, m, Xbounds, prior = "adaptive", r0 = 10)
  expect_equal(model_adaptive$prior, "adaptive")
  expect_equal(model_adaptive$r0, 10)
})

test_that("Basic Functionality: fit.BKP works with different 'kernel' types", {
  model_gaussian <- fit.BKP(X, y, m, Xbounds, kernel = "gaussian")
  expect_equal(model_gaussian$kernel, "gaussian")

  model_matern52 <- fit.BKP(X, y, m, Xbounds, kernel = "matern52")
  expect_equal(model_matern52$kernel, "matern52")
})

test_that("Basic Functionality: fit.BKP works with different 'loss' functions", {
  model_brier <- fit.BKP(X, y, m, Xbounds, loss = "brier")
  expect_equal(model_brier$loss, "brier")

  model_logloss <- fit.BKP(X, y, m, Xbounds, loss = "log_loss")
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
  model <- fit.BKP(X, y, m, Xbounds)
  expect_equal(model$prior, "noninformative")
})

test_that("Default Values: Default 'r0' and 'p0' are used when prior is 'fixed'", {
  # r0=2, p0=0.5 are defaults in function signature
  model <- fit.BKP(X, y, m, Xbounds, prior = "fixed")
  expect_equal(model$r0, 2)
  expect_equal(model$p0, 0.5)
})

test_that("Default Values: Default 'kernel' is 'gaussian'", {
  model <- fit.BKP(X, y, m, Xbounds)
  expect_equal(model$kernel, "gaussian")
})

test_that("Default Values: Default 'loss' is 'brier'", {
  model <- fit.BKP(X, y, m, Xbounds)
  expect_equal(model$loss, "brier")
})
