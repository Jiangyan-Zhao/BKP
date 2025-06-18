test_that("kernel_matrix computes Gaussian kernel correctly", {
  X <- matrix(c(1, 2, 3, 4), nrow = 2)
  theta <- 1
  K <- kernel_matrix(X, X, theta, kernel = "gaussian")
  expect_true(is.matrix(K))
  expect_equal(dim(K), c(2, 2))
  expect_equal(K[1, 1], 1)  # self-distance is zero → exp(0) = 1
  expect_equal(K, t(K), tolerance = 1e-10)  # symmetry
})

test_that("kernel_matrix throws error on mismatched dimensions", {
  X <- matrix(1:6, nrow = 2)
  Xp <- matrix(1:4, nrow = 2)
  expect_error(kernel_matrix(X, Xp, theta = 1), "must have the same number of columns")
})

test_that("kernel_matrix expands scalar and vector theta correctly", {
  X <- matrix(1:6, nrow = 2)
  expect_silent(kernel_matrix(X, theta = 0.5, kernel = "matern32"))
  expect_silent(kernel_matrix(X, theta = c(1, 1), kernel = "matern52"))
  expect_error(kernel_matrix(X, theta = c(1, 2, 3), kernel = "gaussian"),
               "Length of theta must be 1 or equal to number of input dimensions")
})

test_that("kernel_matrix computes values with different kernels and Xprime", {
  X <- matrix(runif(20), ncol = 2)
  Xprime <- matrix(runif(10), ncol = 2)
  theta <- 0.5

  K1 <- kernel_matrix(X, theta = theta, kernel = "gaussian")
  expect_equal(dim(K1), c(nrow(X), nrow(X)))

  K2 <- kernel_matrix(X, Xprime, theta = theta, kernel = "gaussian")
  expect_equal(dim(K2), c(nrow(X), nrow(Xprime)))

  K3 <- kernel_matrix(X, theta = theta, kernel = "matern32")
  expect_true(all(K3 >= 0))

  K4 <- kernel_matrix(X, theta = theta, kernel = "matern52")
  expect_true(all(K4 >= 0))
})

test_that("kernel_matrix throws error on invalid kernel", {
  X <- matrix(1:4, nrow = 2)
  expect_error(kernel_matrix(X, X, theta = 1, kernel = "invalid"), "arg should be one of")
})
