# Test file for the kernel_matrix function

test_that("kernel_matrix handles input validation correctly", {
  # Test with different column numbers for X and Xprime
  expect_error(kernel_matrix(matrix(1:4, 2), matrix(1:6, 2)),
               "'X' and 'Xprime' must have the same number of columns \\(input dimensions\\).")

  # Test with an invalid kernel name
  expect_error(kernel_matrix(matrix(1:4, 2), kernel = "invalid_kernel"),
               "should be one of \"gaussian\", \"matern52\", \"matern32\"", fixed = TRUE)

  # Test with isotropic=TRUE and theta is a vector
  expect_error(kernel_matrix(matrix(1:4, 2), isotropic = TRUE, theta = c(0.5, 0.6)),
               "For isotropic=TRUE, 'theta' must be a scalar.")

  # Test with isotropic=FALSE and theta has wrong length
  expect_error(kernel_matrix(matrix(1:6, 3), isotropic = FALSE, theta = c(0.5, 0.6, 0.7)),
               "For isotropic=FALSE, 'theta' must be scalar or of length equal to ncol\\(X\\).")
})


test_that("kernel_matrix additional validation and branch paths", {
  set.seed(1)
  X <- matrix(runif(6), ncol = 2)
  Xp <- matrix(runif(4), ncol = 2)

  expect_error(kernel_matrix(X = c("a", "b")), "'X' must be numeric or a numeric matrix.")
  X_na <- X
  X_na[1, 1] <- NA
  expect_error(kernel_matrix(X_na), "'X' must contain only finite numeric values.")

  expect_error(kernel_matrix(X, Xprime = c("a", "b")), "'Xprime' must be numeric or a numeric matrix.")
  Xp_na <- Xp
  Xp_na[1, 1] <- NA
  expect_error(kernel_matrix(X, Xp_na), "'Xprime' must contain only finite numeric values.")

  expect_error(kernel_matrix(X, theta = 0), "'theta' must contain only finite numeric and strictly positive values.")
  expect_error(kernel_matrix(X, isotropic = c(TRUE, FALSE)), "'isotropic' must be a single logical value.")

  K_ns <- kernel_matrix(X, Xp, theta = c(0.2, 0.4), kernel = "matern32", isotropic = FALSE)
  expect_equal(dim(K_ns), c(nrow(X), nrow(Xp)))

  K_vec <- kernel_matrix(1:5, theta = 0.5)
  expect_equal(dim(K_vec), c(5, 5))
})


test_that("kernel_matrix returns correct Gaussian values in 1D", {
  X <- matrix(c(0, 1), ncol = 1)
  K <- kernel_matrix(X, theta = 2, kernel = "gaussian")
  expect_equal(K[1, 1], 1)
  expect_equal(K[2, 2], 1)
  expect_equal(K[1, 2], exp(-(1 / 2)^2))
  expect_equal(K[2, 1], K[1, 2])
})

test_that("kernel_matrix rejects non-finite inputs", {
  expect_error(kernel_matrix(matrix(c(0, Inf), ncol = 1)))
  expect_error(kernel_matrix(matrix(c(0, 1), ncol = 1), theta = Inf))
})

test_that("anisotropic scaling is applied dimension-wise", {
  X <- rbind(c(0, 0), c(1, 2))
  K <- kernel_matrix(X, theta = c(1, 2), kernel = "gaussian", isotropic = FALSE)
  expect_equal(K[1, 2], exp(-((1 / 1)^2 + (2 / 2)^2)))
})

test_that("kernel_matrix low-memory engine matches direct Gaussian formula", {
  X <- matrix(seq(0, 1, length.out = 1001), ncol = 1)
  Xp <- matrix(seq(0.25, 0.75, length.out = 1000), ncol = 1)
  theta <- 0.3

  K <- kernel_matrix(X, Xp, theta = theta, kernel = "gaussian")

  idx_i <- c(1, 17, 1001)
  idx_j <- c(1, 503, 1000)
  expected <- exp(-outer(X[idx_i, 1], Xp[idx_j, 1], function(a, b) ((a - b) / theta)^2))
  expect_equal(K[idx_i, idx_j], expected, tolerance = 1e-12)
})

test_that("kernel_matrix loop engine preserves kernel properties", {
  set.seed(1)

  X <- matrix(runif(80), ncol = 4)
  Xprime <- matrix(runif(60), ncol = 4)
  kernels <- c("gaussian", "matern52", "matern32", "wendland")

  for (ker in kernels) {
    K1 <- kernel_matrix(X, theta = 0.3, kernel = ker, isotropic = TRUE)
    K2 <- kernel_matrix(X, Xprime, theta = 0.3, kernel = ker, isotropic = TRUE)

    expect_true(all(is.finite(K1)))
    expect_true(all(is.finite(K2)))
    expect_equal(K1, t(K1), tolerance = 1e-12)
    expect_equal(diag(K1), rep(1, nrow(X)), tolerance = 1e-12)

    K3 <- kernel_matrix(X, theta = c(0.2, 0.3, 0.4, 0.5),
                        kernel = ker, isotropic = FALSE)
    K4 <- kernel_matrix(X, Xprime,
                        theta = c(0.2, 0.3, 0.4, 0.5),
                        kernel = ker, isotropic = FALSE)

    expect_true(all(is.finite(K3)))
    expect_true(all(is.finite(K4)))
    expect_equal(K3, t(K3), tolerance = 1e-12)
    expect_equal(diag(K3), rep(1, nrow(X)), tolerance = 1e-12)
  }
})

test_that("kernel_matrix low-memory engine supports all kernels and lengthscales", {
  set.seed(2)

  X <- matrix(runif(1001 * 2), ncol = 2)
  Xprime <- matrix(runif(1000 * 2), ncol = 2)
  kernels <- c("gaussian", "matern52", "matern32", "wendland")

  idx_i <- c(1, 101, 1001)
  idx_j <- c(1, 500, 1000)

  for (ker in kernels) {
    K_iso <- kernel_matrix(X, Xprime, theta = 0.3, kernel = ker, isotropic = TRUE)
    K_aniso <- kernel_matrix(X, Xprime, theta = c(0.2, 0.5), kernel = ker,
                             isotropic = FALSE)

    expect_true(all(is.finite(K_iso[idx_i, idx_j])))
    expect_true(all(is.finite(K_aniso[idx_i, idx_j])))
    expect_equal(dim(K_iso), c(nrow(X), nrow(Xprime)))
    expect_equal(dim(K_aniso), c(nrow(X), nrow(Xprime)))
  }
})
