# Test file for the fit_BKP function
#
# All other functions are assumed to be available in the package namespace.

# Start of tests
test_that("fit_BKP runs correctly with examples from the documentation", {
  # --- 1D Example ---
  set.seed(123)
  true_pi_fun <- function(x) {
    (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
  }
  n <- 30
  Xbounds <- matrix(c(-2,2), nrow=1)
  X <- lhs(n = n, rect = Xbounds)
  true_pi <- true_pi_fun(X)
  m <- sample(100, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = true_pi)

  expect_no_error({model1 <- fit_BKP(X, y, m, Xbounds=Xbounds)})
  expect_s3_class(model1, "BKP")

  # --- 2D Example ---
  set.seed(123)
  true_pi_fun <- function(X) {
    if(is.null(nrow(X))) X <- matrix(X, nrow=1)
    mn <- 8.6928
    s <- 2.4269
    x1 <- 4*X[,1]- 2
    x2 <- 4*X[,2]- 2
    a <- 1 + (x1 + x2 + 1)^2 *
      (19- 14*x1 + 3*x1^2- 14*x2 + 6*x1*x2 + 3*x2^2)
    b <- 30 + (2*x1- 3*x2)^2 *
      (18- 32*x1 + 12*x1^2 + 48*x2- 36*x1*x2 + 27*x2^2)
    f <- log(a*b)
    f <- (f- mn)/s
    return(pnorm(f))
  }
  n <- 100
  Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
  X <- lhs(n = n, rect = Xbounds)
  true_pi <- true_pi_fun(X)
  m <- sample(100, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = true_pi)

  expect_no_error({model2 <- fit_BKP(X, y, m, Xbounds=Xbounds)})
  expect_s3_class(model2, "BKP")
})

test_that("fit_BKP handles input validation correctly", {
  # Test for missing arguments
  expect_error(fit_BKP(), "Arguments 'X', 'y', and 'm' must be provided.")

  # Define a simple dataset for testing
  n <- 10
  d <- 2
  X_test <- matrix(runif(n * d), nrow = n)
  y_test <- rbinom(n, size = 100, prob = 0.5)
  m_test <- rep(100, n)

  expect_error(fit_BKP(X=X_test, y=y_test), "Arguments 'X', 'y', and 'm' must be provided.")

  # Test for non-numeric or invalid input types
  expect_error(fit_BKP(X=as.character(X_test), y=y_test, m=m_test), "'X' must be a numeric matrix or data frame.")
  expect_error(fit_BKP(X=X_test, y=as.character(y_test), m=m_test), "'y' must be numeric.")

  # Test for dimension mismatches
  expect_error(fit_BKP(X=X_test, y=y_test[1:(n-1)], m=m_test), "'y' must have the same number of rows as 'X'.")

  # Test for invalid y and m values
  expect_error(fit_BKP(X=X_test, y=rep(-1, n), m=m_test), "'y' must be nonnegative.")
  expect_error(fit_BKP(X=X_test, y=y_test, m=rep(0, n)), "'m' must be strictly positive.")
  expect_error(fit_BKP(X=X_test, y=m_test + 1, m=m_test), "Each element of 'y' must be less than or equal to corresponding element of 'm'.")

  # Test for NA values
  X_na <- X_test
  X_na[1, 1] <- NA
  expect_error(fit_BKP(X=X_na, y=y_test, m=m_test), "Missing values are not allowed in 'X', 'y', or 'm'.")

  # Test Xbounds validation
  d <- 2 # Define d for this test block
  expect_error(fit_BKP(X=X_test, y=y_test, m=m_test, Xbounds = 1), "'Xbounds' must be a numeric matrix.")
  expect_error(fit_BKP(X=X_test, y=y_test, m=m_test, Xbounds = matrix(1, nrow = d, ncol = 3)), paste0("'Xbounds' must be a matrix with dimensions d x 2, where d = ", d, "."))
})

test_that("fit_BKP returns a BKP object with correct structure and content", {
  n <- 10
  d <- 2
  X_test <- matrix(runif(n * d), nrow = n)
  y_test <- rbinom(n, size = 100, prob = 0.5)
  m_test <- rep(100, n)

  model <- fit_BKP(X=X_test, y=y_test, m=m_test)

  # Check class and structure
  expect_s3_class(model, "BKP")
  expect_true(is.list(model))
  expect_equal(names(model), c(
    "theta_opt", "kernel", "isotropic", "loss", "loss_min",
    "X", "Xnorm", "Xbounds", "y", "m", "prior", "r0", "p0", "alpha0", "beta0",
    "alpha_n", "beta_n", "ess", "ess_info"
    ))

  # Check content
  expect_equal(model$loss, "brier")
  expect_equal(model$ess, "none")
  expect_equal(model$prior, "noninformative")
  expect_equal(model$X, X_test)
  expect_equal(dim(model$Xnorm), dim(X_test))
})

test_that("fit_BKP uses user-provided theta and skips optimization", {
  user_theta <- 0.1
  n <- 10
  d <- 2
  X_test <- matrix(runif(n * d), nrow = n)
  y_test <- rbinom(n, size = 100, prob = 0.5)
  m_test <- rep(100, n)

  model <- fit_BKP(X=X_test, y=y_test, m=m_test, theta = user_theta)

  expect_equal(model$theta_opt, user_theta)
})


test_that("fit_BKP supports isotropic lengthscale", {
  n <- 10
  d <- 2
  X_test <- matrix(runif(n * d), nrow = n)
  y_test <- rbinom(n, size = 100, prob = 0.5)
  m_test <- rep(100, n)

  model <- fit_BKP(X = X_test, y = y_test, m = m_test, theta = 0.2, isotropic = TRUE)
  expect_equal(model$theta_opt, 0.2)

  expect_error(
    fit_BKP(X = X_test, y = y_test, m = m_test, theta = c(0.2, 0.3), isotropic = TRUE),
    "When isotropic=TRUE, 'theta' must be a scalar."
  )
})

test_that("fit_BKP validates extended argument branches", {
  set.seed(42)
  X <- matrix(runif(18, min = -0.2, max = 1.2), ncol = 3)
  y <- rbinom(nrow(X), size = 5, prob = 0.5)
  m <- rep(5, nrow(X))

  expect_warning(
    expect_s3_class(fit_BKP(X, y, m), "BKP"),
    "Input X does not appear to be normalized to \\[0,1\\]"
  )

  expect_error(
    fit_BKP(X, y, m, Xbounds = matrix(c(0, 0, 1, 0.5, 0.4, 0.4), ncol = 2, byrow = TRUE)),
    "Each row of 'Xbounds' must satisfy lower < upper"
  )

  X <- matrix(runif(18, min = 0, max = 1), ncol = 3)
  y <- rbinom(nrow(X), size = 5, prob = 0.5)
  m <- rep(5, nrow(X))
  expect_error(fit_BKP(X, y, m, r0 = 0), "'r0' must be a positive scalar")
  expect_error(fit_BKP(X, y, m, prior = "fixed", p0 = 1),
               "For fixed prior in BKP, 'p0' must be a scalar in \\(0, 1\\).")
  expect_error(fit_BKP(X, y, m, n_multi_start = 0), "'n_multi_start' must be a positive integer")
  expect_error(fit_BKP(X, y, m, theta = "bad"), "'theta' must be numeric")
  expect_error(fit_BKP(X, y, m, theta = c(0.2, 0.3), isotropic = TRUE), "When isotropic=TRUE")
  expect_error(fit_BKP(X, y, m, theta = c(0.2, 0.3), isotropic = FALSE), "length 3")
  expect_error(fit_BKP(X, y, m, theta = -0.2), "'theta' must be strictly positive")
  expect_error(fit_BKP(X, y, m, isotropic = c(TRUE, FALSE)), "'isotropic' must be a single logical value")
})

test_that("ess = 'none' matches the standard BKP posterior update", {
  X <- matrix(c(0.05, 0.20, 0.45, 0.70, 0.90), ncol = 1)
  m <- c(10, 20, 15, 25, 30)
  y <- c(2, 12, 7, 15, 21)

  model <- fit_BKP(X, y, m, theta = 0.3, ess = "none")
  K <- kernel_matrix(model$Xnorm, theta = model$theta_opt,
                     kernel = model$kernel, isotropic = model$isotropic)
  prior_par <- get_prior(prior = model$prior, model = "BKP",
                         r0 = model$r0, p0 = model$p0,
                         y = model$y, m = model$m, K = K)

  expect_equal(model$ess, "none")
  expect_equal(model$alpha_n, prior_par$alpha0 + as.vector(K %*% model$y))
  expect_equal(model$beta_n, prior_par$beta0 + as.vector(K %*% (model$m - model$y)))

  pred_train <- predict(model)
  pred_new <- predict(model, Xnew = matrix(c(0.15, 0.55), ncol = 1))
  expect_equal(pred_train$ess, "none")
  expect_equal(pred_new$ess, "none")
  expect_equal(pred_train$ess_info$scale, rep(1, length(m)))
  expect_equal(pred_new$ess_info$scale, rep(1, 2))
})

test_that("Shepard ESS targets observed trial sizes at training points", {
  X <- matrix(c(0.02, 0.18, 0.44, 0.73, 0.96), ncol = 1)
  m <- c(8, 17, 11, 23, 31)
  y <- c(3, 8, 4, 15, 22)

  model <- fit_BKP(X, y, m, theta = 0.25, ess = "shepard")
  data_size <- model$alpha_n + model$beta_n - model$alpha0 - model$beta0

  expect_equal(as.vector(data_size), as.numeric(m), tolerance = 1e-10)
  expect_equal(model$ess_info$m_shepard, as.numeric(m), tolerance = 1e-10)
})

test_that("Shepard ESS preserves the kernel-weighted empirical success proportion", {
  X <- matrix(c(0.10, 0.15, 0.35, 0.60, 0.85, 0.95), ncol = 2, byrow = TRUE)
  m <- c(10, 12, 18)
  y <- c(4, 9, 7)

  model <- fit_BKP(X, y, m, theta = 0.4, ess = "shepard")
  K <- kernel_matrix(model$Xnorm, theta = model$theta_opt,
                     kernel = model$kernel, isotropic = model$isotropic)
  raw_success <- as.vector(K %*% model$y)
  raw_trials <- as.vector(K %*% model$m)
  c_scale <- model$ess_info$scale

  positive <- raw_trials > 0
  expect_equal(
    raw_success[positive] / raw_trials[positive],
    (c_scale[positive] * raw_success[positive]) /
      (c_scale[positive] * raw_trials[positive]),
    tolerance = 1e-12
  )
})

test_that("predict works for Shepard ESS at new points without NA values", {
  X <- matrix(c(0.05, 0.25, 0.50, 0.80), ncol = 1)
  m <- c(10, 20, 15, 30)
  y <- c(2, 11, 8, 21)
  Xnew <- matrix(c(0.10, 0.40, 0.65, 0.90), ncol = 1)

  model <- fit_BKP(X, y, m, theta = 0.35, ess = "shepard")
  pred <- predict(model, Xnew = Xnew)

  expect_false(anyNA(pred$alpha_n))
  expect_false(anyNA(pred$beta_n))
  expect_false(anyNA(pred$mean))
  expect_true(all(is.finite(pred$alpha_n)))
  expect_true(all(is.finite(pred$beta_n)))
  expect_true(all(is.finite(pred$mean)))
})

test_that("Shepard ESS rejects duplicated input locations", {
  X <- matrix(c(0.10, 0.20, 0.10, 0.20, 0.75, 0.90), ncol = 2, byrow = TRUE)
  m <- c(10, 20, 15)
  y <- c(4, 9, 7)

  expect_error(
    fit_BKP(X, y, m, theta = 0.3, ess = "shepard"),
    "requires unique input locations"
  )
})
