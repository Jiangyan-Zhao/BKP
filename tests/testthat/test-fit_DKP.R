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
  expect_equal(names(model), c(
    "theta_opt", "kernel", "isotropic", "loss", "loss_min",
    "X", "Xnorm", "Xbounds", "Y", "prior", "r0",
    "p0", "alpha0", "alpha_n", "ess", "ess_info"))

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
    fit_DKP(X=X_test, Y=Y_test, theta = user_theta, isotropic = FALSE),
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


test_that("fit_DKP supports isotropic lengthscale", {
  n <- 12
  d <- 2
  X_test <- matrix(runif(n * d), nrow = n)
  Y_test <- t(sapply(1:n, function(i) rmultinom(1, size = 50, prob = c(0.2, 0.3, 0.5))))

  model <- fit_DKP(X = X_test, Y = Y_test, theta = 0.2, isotropic = TRUE)
  expect_equal(model$theta_opt, 0.2)

  expect_error(
    fit_DKP(X = X_test, Y = Y_test, theta = c(0.2, 0.3), isotropic = TRUE),
    "When isotropic=TRUE, 'theta' must be a scalar."
  )
})

test_that("fit_DKP validates extended argument branches", {
  set.seed(43)
  X <- matrix(runif(18, min = -0.2, max = 1.2), ncol = 3)
  probs <- c(0.2, 0.3, 0.5)
  Y <- t(sapply(seq_len(nrow(X)), function(i) rmultinom(1, size = 8, prob = probs)))

  expect_warning(
    expect_s3_class(fit_DKP(X, Y), "DKP"),
    "Input X does not appear to be normalized to \\[0,1\\]"
  )

  expect_error(
    fit_DKP(X, Y, Xbounds = matrix(c(0, 0, 1, 0.5, 0.4, 0.4), ncol = 2, byrow = TRUE)),
    "Each row of 'Xbounds' must satisfy lower < upper"
  )

  X <- matrix(runif(18, min = 0, max = 1), ncol = 3)
  probs <- c(0.2, 0.3, 0.5)
  Y <- t(sapply(seq_len(nrow(X)), function(i) rmultinom(1, size = 8, prob = probs)))
  expect_error(fit_DKP(X, Y, r0 = 0), "'r0' must be a positive scalar")
  expect_error(fit_DKP(X, Y, prior = "fixed", p0 = c(0.7, -0.2, 0.5)),
               "For fixed prior in DKP, 'p0' must be a nonnegative finite numeric vector of length equal to the number of classes and sum to 1.")
  expect_error(fit_DKP(X, Y, prior = "fixed", p0 = c(0.5, 0.5)),
               "For fixed prior in DKP, 'p0' must be a nonnegative finite numeric vector of length equal to the number of classes and sum to 1.")
  expect_error(fit_DKP(X, Y, n_multi_start = 0), "'n_multi_start' must be a positive integer")
  expect_error(fit_DKP(X, Y, theta = "bad"), "'theta' must be numeric")
  expect_error(fit_DKP(X, Y, theta = c(0.2, 0.3), isotropic = TRUE), "When isotropic=TRUE")
  expect_error(fit_DKP(X, Y, theta = c(0.2, 0.3), isotropic = FALSE), "length 3")
  expect_error(fit_DKP(X, Y, theta = -0.2), "'theta' must be strictly positive")
  expect_error(fit_DKP(X, Y, isotropic = c(TRUE, FALSE)), "'isotropic' must be a single logical value")
})


test_that("DKP ess = 'none' matches the standard posterior update", {
  X <- matrix(c(0.05, 0.35, 0.70, 0.95), ncol = 1)
  Y <- matrix(c(3, 1, 2,
                1, 4, 2,
                2, 2, 5,
                5, 1, 1), ncol = 3, byrow = TRUE)

  model <- fit_DKP(X, Y, theta = 0.3, prior = "fixed", r0 = 2,
                   p0 = c(1/3, 1/3, 1/3), ess = "none")
  K <- kernel_matrix(model$Xnorm, theta = model$theta_opt,
                     kernel = model$kernel, isotropic = model$isotropic)

  expect_equal(model$ess, "none")
  expect_equal(model$alpha_n, model$alpha0 + as.matrix(K %*% Y), tolerance = 0)
})

test_that("DKP Shepard ESS calibrates training-point total data contribution", {
  X <- matrix(c(0.05, 0.35, 0.70, 0.95), ncol = 1)
  Y <- matrix(c(3, 1, 2,
                1, 4, 2,
                2, 2, 5,
                5, 1, 1), ncol = 3, byrow = TRUE)
  m <- rowSums(Y)

  model <- fit_DKP(X, Y, theta = 0.25, prior = "fixed", r0 = 2,
                   p0 = c(1/3, 1/3, 1/3), ess = "shepard")

  expect_equal(rowSums(model$alpha_n - model$alpha0), as.numeric(m), tolerance = 1e-8)
  expect_equal(model$ess_info$m_shepard, as.numeric(m), tolerance = 1e-10)
})

test_that("DKP Shepard ESS preserves kernel-weighted empirical class proportions", {
  X <- matrix(c(0.05, 0.35, 0.70, 0.95), ncol = 1)
  Y <- matrix(c(3, 1, 2,
                1, 4, 2,
                2, 2, 5,
                5, 1, 1), ncol = 3, byrow = TRUE)

  model <- fit_DKP(X, Y, theta = 0.4, prior = "fixed", r0 = 2,
                   p0 = c(1/3, 1/3, 1/3), ess = "shepard")
  K <- kernel_matrix(model$Xnorm, theta = model$theta_opt,
                     kernel = model$kernel, isotropic = model$isotropic)
  raw_counts <- as.matrix(K %*% Y)
  raw_totals <- rowSums(raw_counts)
  scaled_counts <- sweep(raw_counts, 1L, model$ess_info$scale, "*")
  scaled_totals <- rowSums(scaled_counts)
  positive <- raw_totals > 0

  expect_equal(scaled_counts[positive, , drop = FALSE] / scaled_totals[positive],
               raw_counts[positive, , drop = FALSE] / raw_totals[positive],
               tolerance = 1e-10)
})

test_that("DKP Shepard ESS prediction on new points has no NA values", {
  X <- matrix(c(0.05, 0.35, 0.70, 0.95), ncol = 1)
  Y <- matrix(c(3, 1, 2,
                1, 4, 2,
                2, 2, 5,
                5, 1, 1), ncol = 3, byrow = TRUE)

  model <- fit_DKP(X, Y, theta = 0.35, ess = "shepard")
  pred <- predict(model, Xnew = matrix(c(0.2, 0.6, 0.9), ncol = 1))

  expect_false(anyNA(pred$alpha_n))
  expect_false(anyNA(pred$mean))
  expect_equal(pred$ess, "shepard")
})

test_that("DKP Shepard ESS rejects duplicated input locations", {
  X <- matrix(c(0.1, 0.1, 0.7), ncol = 1)
  Y <- matrix(c(2, 1, 1,
                1, 3, 1,
                2, 2, 4), ncol = 3, byrow = TRUE)

  expect_error(
    fit_DKP(X, Y, theta = 0.3, ess = "shepard"),
    "requires unique input locations"
  )
})

test_that("DKP Shepard ESS works with optimized and fixed theta", {
  set.seed(12)
  X <- matrix(seq(0.05, 0.95, length.out = 6), ncol = 1)
  Y <- matrix(c(3, 1, 2,
                1, 4, 2,
                2, 2, 5,
                5, 1, 1,
                3, 3, 2,
                1, 2, 6), ncol = 3, byrow = TRUE)

  optimized <- fit_DKP(X, Y, theta = NULL, n_multi_start = 2, ess = "shepard")
  fixed <- fit_DKP(X, Y, theta = 0.3, ess = "shepard")

  expect_s3_class(optimized, "DKP")
  expect_s3_class(fixed, "DKP")
  expect_equal(optimized$ess, "shepard")
  expect_equal(fixed$ess, "shepard")
})
