# Test file for the fit_DKP function
#
# All other functions are assumed to be available in the package namespace.

# Start of tests
test_that("fit_DKP runs correctly with examples from the documentation", {
  # --- 1D Example ---
  set.seed(123)

  # Define true class probability function (3-class)
  true_pi_fun <- function(X) {
    p1 <- 1/(1+exp(-3*X))
    p2 <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
    return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
  }

  n <- 30
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- lhs(n = n, rect = Xbounds)
  true_pi <- true_pi_fun(X)
  m <- sample(150, n, replace = TRUE)
  Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))

  # Fit DKP model
  expect_no_error({model1 <- fit_DKP(X, Y, Xbounds = Xbounds)})
  expect_s3_class(model1, "DKP")

  # --- 2D Example ---
  set.seed(123)

  # Define latent function and transform to 3-class probabilities
  true_pi_fun <- function(X) {
    if (is.null(nrow(X))) X <- matrix(X, nrow = 1)
    mn <- 8.6928; s <- 2.4269
    x1 <- 4 * X[,1] - 2
    x2 <- 4 * X[,2] - 2
    a <- 1 + (x1 + x2 + 1)^2 *
      (19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2)
    b <- 30 + (2*x1 - 3*x2)^2 *
      (18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2)
    f <- (log(a*b)- mn)/s
    p1 <- pnorm(f) # Transform to probability
    p2 <- sin(pi * X[,1]) * sin(pi * X[,2])
    return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
  }

  n <- 100
  Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
  X <- lhs(n = n, rect = Xbounds)
  true_pi <- true_pi_fun(X)
  m <- sample(150, n, replace = TRUE)
  Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))

  # Fit DKP model
  expect_no_error({model2 <- fit_DKP(X, Y, Xbounds = Xbounds)})
  expect_s3_class(model2, "DKP")
})

test_that("fit_DKP handles input validation correctly", {
  # Test for missing arguments
  expect_error(fit_DKP(), "Arguments 'X' and 'Y' must be provided.")

  # Define a simple dataset for testing
  n <- 10
  d <- 2
  X_test <- matrix(runif(n * d), nrow = n)
  Y_test <- t(sapply(1:n, function(i) rmultinom(1, size = 100, prob = c(0.5, 0.5))))

  # Test for invalid input types
  expect_error(fit_DKP(X=as.character(X_test), Y=Y_test), "'X' must be a numeric matrix or data frame.")
  expect_error(fit_DKP(X=X_test, Y=as.character(Y_test)), "'Y' must be a numeric matrix or data frame.")

  # Test for non-numeric input
  expect_error(fit_DKP(X=matrix("a", nrow=n, ncol=d), Y=Y_test), "'X' must contain numeric values only.")
  expect_error(fit_DKP(X=X_test, Y=matrix("a", nrow=n, ncol=2)), "'Y' must contain numeric values only.")

  # Test for dimension mismatches
  expect_error(fit_DKP(X=X_test, Y=Y_test[1:(n-1),]), "Number of rows in 'Y' must match number of rows in 'X'.")

  # Test for invalid Y values
  Y_neg <- Y_test
  Y_neg[1, 1] <- -1
  expect_error(fit_DKP(X=X_test, Y=Y_neg), "'Y' must be nonnegative counts or frequencies.")

  # Test for NA values
  X_na <- X_test
  X_na[1, 1] <- NA
  expect_error(fit_DKP(X=X_na, Y=Y_test), "Missing values are not allowed in 'X' or 'Y'.")

  # Test Xbounds validation
  d <- 2 # Define d for this test block
  expect_warning(expect_error(fit_DKP(X=X_test, Y=Y_test, Xbounds = 1), "'Xbounds' must be a numeric matrix."))
  expect_warning(expect_error(fit_DKP(X=X_test, Y=Y_test, Xbounds = matrix(1, nrow = d, ncol = 3)), paste0("'Xbounds' must be a matrix with dimensions d x 2, where d = ", d, ".")))
})

test_that("fit_DKP returns a DKP object with correct structure and content", {
  n <- 10
  d <- 2
  X_test <- matrix(runif(n * d), nrow = n)
  Y_test <- t(sapply(1:n, function(i) rmultinom(1, size = 100, prob = c(0.5, 0.5))))

  # Use withCallingHandlers to capture the warning and get the return value
  model_warning <- NULL
  model <- withCallingHandlers(
    fit_DKP(X=X_test, Y=Y_test),
    warning = function(w) {
      model_warning <<- w
      invokeRestart("muffleWarning")
    }
  )

  # Expect the specific warning and check model class
  expect_s3_class(model, "DKP")
  expect_s3_class(model_warning, "warning")
  expect_equal(model_warning$message, "For binary data, consider using the BKP model instead of DKP.")

  # Check class and structure
  expect_true(is.list(model))
  expect_equal(names(model), c("theta_opt", "kernel", "loss", "loss_min", "X", "Xnorm", "Xbounds", "Y", "prior", "r0", "p0", "alpha0", "alpha_n"))

  # Check content
  expect_equal(model$loss, "brier")
  expect_equal(model$prior, "noninformative")
  expect_equal(model$X, X_test)
  expect_equal(dim(model$Xnorm), dim(X_test))
})

test_that("fit_DKP uses user-provided theta and skips optimization", {
  user_theta <- c(0.1, 0.5)
  n <- 10
  d <- 2
  X_test <- matrix(runif(n * d), nrow = n)
  Y_test <- t(sapply(1:n, function(i) rmultinom(1, size = 100, prob = c(0.5, 0.5))))

  # Use withCallingHandlers to capture the warning and get the return value
  model_warning <- NULL
  model <- withCallingHandlers(
    fit_DKP(X=X_test, Y=Y_test, theta = user_theta),
    warning = function(w) {
      model_warning <<- w
      invokeRestart("muffleWarning")
    }
  )

  # Expect the specific warning
  expect_s3_class(model_warning, "warning")
  expect_equal(model_warning$message, "For binary data, consider using the BKP model instead of DKP.")

  expect_equal(model$theta_opt, user_theta)
})
