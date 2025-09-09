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
  expect_equal(names(model), c("theta_opt", "kernel", "loss", "loss_min", "X", "Xnorm", "Xbounds", "y", "m", "prior", "r0", "p0", "alpha0", "beta0", "alpha_n", "beta_n"))

  # Check content
  expect_equal(model$loss, "brier")
  expect_equal(model$prior, "noninformative")
  expect_equal(model$X, X_test)
  expect_equal(dim(model$Xnorm), dim(X_test))
})

test_that("fit_BKP uses user-provided theta and skips optimization", {
  user_theta <- c(0.1, 0.5)
  n <- 10
  d <- 2
  X_test <- matrix(runif(n * d), nrow = n)
  y_test <- rbinom(n, size = 100, prob = 0.5)
  m_test <- rep(100, n)

  model <- fit_BKP(X=X_test, y=y_test, m=m_test, theta = user_theta)

  expect_equal(model$theta_opt, user_theta)
})
