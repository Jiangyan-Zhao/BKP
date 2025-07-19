
# -------------------------- Test Setup ---------------------------
# Define some common test data
X_test <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE) # 4x2 matrix
Xprime_test <- matrix(c(0.5, 0.5, 0, 0), ncol = 2, byrow = TRUE) # 2x2 matrix

# -------------------------- Test Context: Basic Functionality ---------------------------

test_that("Basic Functionality: kernel_matrix returns a numeric matrix", {
  K <- kernel_matrix(X_test)
  expect_true(is.matrix(K))
  expect_type(K, "double")
})

test_that("default Xprime=X results in symmetric matrix with correct dimensions", {
  K <- kernel_matrix(X_test)
  expect_equal(dim(K), c(nrow(X_test), nrow(X_test)))
  expect_true(isSymmetric(K))
})

test_that("Xprime different from X results in correct dimensions", {
  K <- kernel_matrix(X_test, Xprime_test)
  expect_equal(dim(K), c(nrow(X_test), nrow(Xprime_test)))
  expect_false(isSymmetric(K)) # Should not be symmetric generally
})

test_that("gaussian kernel works correctly and produces values between 0 and 1", {
  K <- kernel_matrix(X_test, kernel = "gaussian", theta = 1)
  expect_true(all(K >= 0 & K <= 1))
  # Expect diagonal to be 1 for Gaussian kernel when Xprime is X
  expect_equal(as.numeric(diag(K)), rep(1, nrow(X_test)))
})

test_that("matern52 kernel works correctly and produces non-negative values", {
  K <- kernel_matrix(X_test, kernel = "matern52", theta = 1)
  expect_true(all(K >= 0))
  expect_equal(as.numeric(diag(K)), rep(1, nrow(X_test)))
})

test_that("matern32 kernel works correctly and produces non-negative values", {
  K <- kernel_matrix(X_test, kernel = "matern32", theta = 1)
  expect_true(all(K >= 0))
  expect_equal(as.numeric(diag(K)), rep(1, nrow(X_test)))
})

test_that("anisotropic=TRUE applies different lengthscales per dimension", {
  theta_aniso <- c(0.1, 10)
  X_small <- matrix(c(0,0, 1,0, 0,1), ncol=2, byrow=TRUE) # (0,0), (1,0), (0,1)
  # Distance between (0,0) and (1,0) should be large for theta[1]=0.1
  # Distance between (0,0) and (0,1) should be small for theta[2]=10
  K_aniso <- kernel_matrix(X_small, theta = theta_aniso, kernel = "gaussian", anisotropic = TRUE)

  # Manually check a few values based on expected behavior
  # (0,0) vs (1,0) -> scaled dist for x1: (1-0)/0.1 = 10, for x2: (0-0)/10 = 0 -> dist_sq = 100
  # (0,0) vs (0,1) -> scaled dist for x1: (0-0)/0.1 = 0, for x2: (1-0)/10 = 0.1 -> dist_sq = 0.01
  expect_equal(K_aniso[1,2], exp(-100)) # Similarity for (0,0) and (1,0)
  expect_equal(K_aniso[1,3], exp(-0.01)) # Similarity for (0,0) and (0,1)
})

test_that("anisotropic=FALSE uses a single scalar lengthscale for all dimensions", {
  theta_iso <- 1
  X_small <- matrix(c(0,0, 1,0, 0,1), ncol=2, byrow=TRUE)
  K_iso <- kernel_matrix(X_small, theta = theta_iso, kernel = "gaussian", anisotropic = FALSE)

  # Distances:
  # (0,0) to (1,0) -> sqrt(1^2 + 0^2) = 1 -> K = exp(-1^2)
  # (0,0) to (0,1) -> sqrt(0^2 + 1^2) = 1 -> K = exp(-1^2)
  expect_equal(K_iso[1,2], exp(-1))
  expect_equal(K_iso[1,3], exp(-1))
})

# -------------------------- Test Context: Input Validation ---------------------------

test_that("Input Validation: X and Xprime must have same number of columns", {
  X_mismatched <- matrix(runif(10), ncol = 1)
  expect_error(kernel_matrix(X_test, X_mismatched))
})

test_that("Input Validation: theta length must be 1 or equal to d when anisotropic=TRUE", {
  X_d3 <- matrix(runif(30), ncol=3)
  expect_error(kernel_matrix(X_d3, theta = c(0.1, 0.2), anisotropic = TRUE)) # length != d
  expect_silent(kernel_matrix(X_d3, theta = 0.1, anisotropic = TRUE)) # scalar is ok
  expect_silent(kernel_matrix(X_d3, theta = c(0.1, 0.2, 0.3), anisotropic = TRUE)) # length == d is ok
})

test_that("Input Validation: theta must be a scalar when anisotropic=FALSE", {
  expect_error(kernel_matrix(X_test, theta = c(0.1, 0.2), anisotropic = FALSE)) # vector not ok
  expect_silent(kernel_matrix(X_test, theta = 0.1, anisotropic = FALSE)) # scalar is ok
})

test_that("Input Validation: unsupported kernel type throws error", {
  expect_error(kernel_matrix(X_test, kernel = "invalid_kernel"))
})

test_that("Input Validation: vector inputs X and Xprime are handled", {
  X_vec <- c(1,2,3)
  Xprime_vec <- c(1.5, 2.5)
  K_vec <- kernel_matrix(X_vec, Xprime_vec, kernel = "gaussian")
  expect_true(is.matrix(K_vec))
  expect_equal(dim(K_vec), c(length(X_vec), length(Xprime_vec)))
})
