# -------------------------- Test Setup ---------------------------
set.seed(123)
n_test <- 50
d_test <- 2
q_test <- 3 # Number of classes
X_test <- matrix(runif(n_test * d_test, -5, 5), n_test, d_test)
# Generate multinomial counts, ensuring sum for each row
Y_test_m <- sample(50, n_test, replace = TRUE) # Total trials for DKP
Y_test <- t(sapply(1:n_test, function(i) rmultinom(1, size = Y_test_m[i], prob = rep(1/q_test, q_test))))

# -------------------------- Test Context: Basic Functionality ---------------------------

test_that("Basic Functionality: fit.DKP returns a 'DKP' class object with expected elements", {
  model <- fit.DKP(X_test, Y_test)

  expect_s3_class(model, "DKP")
  expect_type(model, "list")

  # 检查返回列表中的关键元素
  expect_named(model, c("theta_opt", "kernel", "loss", "loss_min",
                        "X", "Xnorm", "Xbounds", "Y",
                        "prior", "r0", "p0",
                        "alpha0", "alpha_n"))

  expect_true(all(model$X == X_test))
  expect_true(all(model$Y == Y_test))
  expect_equal(length(model$theta_opt), d_test)
  expect_type(model$loss_min, "double")
  expect_equal(ncol(model$Y), q_test)
  expect_equal(ncol(model$alpha0), q_test)
  expect_equal(ncol(model$alpha_n), q_test)
})

test_that("Basic Functionality: fit.DKP works with different 'prior' types", {
  model_noninformative <- fit.DKP(X_test, Y_test, prior = "noninformative")
  expect_equal(model_noninformative$prior, "noninformative")

  # For 'fixed' prior, p0 is required for DKP
  p0_fixed <- rep(1/q_test, q_test)
  model_fixed <- fit.DKP(X_test, Y_test, prior = "fixed", r0 = 5, p0 = p0_fixed)
  expect_equal(model_fixed$prior, "fixed")
  expect_equal(model_fixed$r0, 5)
  expect_equal(model_fixed$p0, p0_fixed)

  model_adaptive <- fit.DKP(X_test, Y_test, prior = "adaptive", r0 = 10)
  expect_equal(model_adaptive$prior, "adaptive")
  expect_equal(model_adaptive$r0, 10)
})

test_that("Basic Functionality: fit.DKP works with different 'kernel' types", {
  model_gaussian <- fit.DKP(X_test, Y_test, kernel = "gaussian")
  expect_equal(model_gaussian$kernel, "gaussian")

  model_matern52 <- fit.DKP(X_test, Y_test, kernel = "matern52")
  expect_equal(model_matern52$kernel, "matern52")
})

test_that("Basic Functionality: fit.DKP works with different 'loss' functions", {
  model_brier <- fit.DKP(X_test, Y_test, loss = "brier")
  expect_equal(model_brier$loss, "brier")

  model_logloss <- fit.DKP(X_test, Y_test, loss = "log_loss")
  expect_equal(model_logloss$loss, "log_loss")
})

# -------------------------- Test Context: Input Validation ---------------------------

test_that("Input Validation: 'Y' must match nrow(X)", {
  Y_invalid <- Y_test[1:10, ]
  expect_error(fit.DKP(X_test, Y_invalid))
})

test_that("Input Validation: 'Y' must be non-negative", {
  Y_invalid <- Y_test
  Y_invalid[1, 1] <- -1
  expect_error(fit.DKP(X_test, Y_invalid))
})

test_that("Input Validation: 'Xbounds' must be a d x 2 matrix", {
  expect_error(fit.DKP(X_test, Y_test, Xbounds = matrix(c(0,1), ncol=1))) # Wrong dimensions
  expect_error(fit.DKP(X_test, Y_test, Xbounds = matrix(c(0,1,0,1,0,1), ncol=3))) # Wrong dimensions
})

test_that("Input Validation: p0 is required when prior is 'fixed'", {
  expect_error(fit.DKP(X_test, Y_test, prior = "fixed", r0 = 5, p0 = NULL)) # p0 is NULL by default
  expect_error(fit.DKP(X_test, Y_test, prior = "fixed", p0 = c(0.1))) # p0 wrong length
})

# -------------------------- Test Context: Default Values ---------------------------

test_that("Default Values: Default 'prior' is 'noninformative'", {
  model <- fit.DKP(X_test, Y_test)
  expect_equal(model$prior, "noninformative")
})

test_that("Default Values: Default 'r0' is used when prior is 'fixed' or 'adaptive'", {
  model_fixed <- fit.DKP(X_test, Y_test, prior = "fixed", p0 = rep(1/q_test, q_test))
  expect_equal(model_fixed$r0, 2) # r0=2 is default

  model_adaptive <- fit.DKP(X_test, Y_test, prior = "adaptive")
  expect_equal(model_adaptive$r0, 2) # r0=2 is default
})

test_that("Default Values: Default 'kernel' is 'gaussian'", {
  model <- fit.DKP(X_test, Y_test)
  expect_equal(model$kernel, "gaussian")
})

test_that("Default Values: Default 'loss' is 'brier'", {
  model <- fit.DKP(X_test, Y_test)
  expect_equal(model$loss, "brier")
})

# -------------------------- Test Context: Edge Cases ---------------------------

test_that("Edge Cases: fit.DKP handles 1D input (d=1)", {
  X_1d <- matrix(runif(n_test, 0, 1), ncol = 1)
  model <- fit.DKP(X_1d, Y_test)
  expect_s3_class(model, "DKP")
  expect_equal(ncol(model$X), 1)
  expect_equal(length(model$theta_opt), 1)
})

test_that("Edge Cases: fit.DKP handles small number of observations (n)", {
  X_small <- matrix(runif(5 * d_test), 5, d_test)
  Y_small <- t(sapply(1:5, function(i) rmultinom(1, size = 10, prob = rep(1/q_test, q_test))))
  model <- fit.DKP(X_small, Y_small)
  expect_s3_class(model, "DKP")
  expect_equal(nrow(model$X), 5)
})
