test_that("kernel_matrix computes Gaussian kernel correctly", {
  X <- matrix(c(1, 2, 3, 4), nrow = 2)
  theta <- 1
  K <- kernel_matrix(X, X, theta, kernel = "gaussian")
  expect_true(is.matrix(K))
  expect_equal(dim(K), c(2, 2))
  expect_equal(K[1, 1], 1)  # self-distance is zero → exp(0) = 1
  expect_equal(K, t(K), tolerance = 1e-10)  # symmetry
})
