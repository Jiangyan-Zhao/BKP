# Test file for the fit_BKP function
#
# All other functions are assumed to be available in the package namespace.

# Start of tests
test_that("fit_BKP runs correctly with examples from the documentation", {
  # --- 1D Example ---
  set.seed(123)

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

  # Fit BKP model
  expect_no_error({model1 <- fit_BKP(X, y, m, Xbounds = Xbounds)})
  expect_s3_class(model1, "BKP")
  expect_equal(dim(model1$Xnorm), dim(X))

  # --- 2D Example ---
  set.seed(123)

  # Define 2D latent function and probability transformation
  true_pi_fun <- function(X) {
    if (is.null(nrow(X))) X <- matrix(X, nrow = 1)
    m_val <- 8.6928
    s <- 2.4269
    x1 <- 4 * X[,1] - 2
    x2 <- 4 * X[,2] - 2
    a <- 1 + (x1 + x2 + 1)^2 *
      (19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2)
    b <- 30 + (2*x1 - 3*x2)^2 *
      (18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2)
    f <- log(a * b)
    f <- (f - m_val) / s
    return(pnorm(f))
  }

  n <- 100
  Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
  X <- tgp::lhs(n = n, rect = Xbounds)
  true_pi <- true_pi_fun(X)
  m <- sample(100, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = true_pi)

  # Fit BKP model
  expect_no_error({model2 <- fit_BKP(X, y, m, Xbounds = Xbounds)})
  expect_s3_class(model2, "BKP")
  expect_equal(dim(model2$Xnorm), dim(X))
})

test_that("fit_BKP handles input validation correctly", {
  # Define a simple dataset for testing
  n <- 10
  d <- 2
  X_test <- matrix(runif(n * d), nrow = n)
  m_test <- rep(100, n)
  y_test <- rbinom(n, size = m_test, prob = 0.5)

  # Test for missing arguments
  expect_error(fit_BKP(), "Arguments 'X', 'y', and 'm' must be provided.")
  expect_error(fit_BKP(X = X_test), "Arguments 'X', 'y', and 'm' must be provided.")

  # Test for invalid input types
  expect_error(fit_BKP(X = as.character(X_test), y = y_test, m = m_test), "'X' must be a numeric matrix or data frame.")
  expect_error(fit_BKP(X = X_test, y = as.character(y_test), m = m_test), "'y' must be numeric.")
  expect_error(fit_BKP(X = X_test, y = y_test, m = as.character(m_test)), "'m' must be numeric.")

  # Test for dimension mismatches
  expect_error(fit_BKP(X = X_test, y = y_test[1:(n-1)], m = m_test), "'y' must have the same number of rows as 'X'.")

  # Test for invalid y and m values
  y_neg <- y_test
  y_neg[1] <- -1
  expect_error(fit_BKP(X = X_test, y = y_neg, m = m_test), "'y' must be nonnegative.")

  m_zero <- m_test
  m_zero[1] <- 0
  expect_error(fit_BKP(X = X_test, y = y_test, m = m_zero), "'m' must be strictly positive.")

  y_gt_m <- y_test
  y_gt_m[1] <- m_test[1] + 1
  expect_error(fit_BKP(X = X_test, y = y_gt_m, m = m_test), "Each element of 'y' must be less than or equal to corresponding element of 'm'.")

  # Test for NA values
  X_na <- X_test
  X_na[1, 1] <- NA
  expect_error(fit_BKP(X = X_na, y = y_test, m = m_test), "Missing values are not allowed in 'X', 'y', or 'm'.")

  # Test Xbounds validation
  expect_error(fit_BKP(X = X_test, y = y_test, m = m_test, Xbounds = 1), "'Xbounds' must be a numeric matrix.")
  expect_error(fit_BKP(X = X_test, y = y_test, m = m_test, Xbounds = matrix(1, nrow = d, ncol = 3)), "'Xbounds' must be a matrix with dimensions d x 2, where d = 2.")
})

test_that("fit_BKP returns a BKP object with correct structure and content", {
  n <- 10
  d <- 2
  X_test <- matrix(runif(n * d), nrow = n)
  m_test <- rep(100, n)
  y_test <- rbinom(n, size = m_test, prob = 0.5)

  model <- fit_BKP(X = X_test, y = y_test, m = m_test)

  # Check class and structure
  expect_s3_class(model, "BKP")
  expect_true(is.list(model))
  expect_equal(names(model), c("theta_opt", "kernel", "loss", "loss_min", "X", "Xnorm", "Xbounds", "y", "m", "prior", "r0", "p0", "alpha0", "beta0", "alpha_n", "beta_n"))

  # Check content
  expect_equal(model$loss, "brier")
  expect_equal(model$prior, "noninformative")
  expect_equal(model$X, X_test)
  expect_equal(model$y, as.matrix(y_test))
  expect_equal(model$m, as.matrix(m_test))
  expect_equal(dim(model$Xnorm), dim(X_test))
})

test_that("fit_BKP uses user-provided theta and skips optimization", {
  user_theta <- c(0.1, 0.5)
  n <- 10
  d <- 2
  X_test <- matrix(runif(n * d), nrow = n)
  m_test <- rep(100, n)
  y_test <- rbinom(n, size = m_test, prob = 0.5)

  model <- fit_BKP(X = X_test, y = y_test, m = m_test, theta = user_theta)

  expect_equal(model$theta_opt, user_theta)
})
