# This test file validates the helper functions in utils.R, which are used for plotting.
# We check if the functions can be called with valid data without throwing an error.

# Load the testthat library for unit testing.
library(testthat)

# The plotting functions rely on these packages.
# They should be loaded for the tests to run.
library(lattice)
library(gridExtra)

# Create mock data that simulates the input for the plotting functions.
# This data mimics the structure of what the main plot methods would pass.
set.seed(123)
n <- 50
d <- 2
q <- 3

# Mock input data (X) and responses (y, Y)
X <- matrix(runif(n * d, -1, 1), ncol = d)
y_bin <- rbinom(n, 1, 0.5)
Y_multi <- t(sapply(1:n, function(i) rmultinom(1, size = 1, prob = c(0.2, 0.5, 0.3))))

# Mock grid data, simulating the output from a 'predict' function.
# This data contains columns for the plot's variables (e.g., Mean, Variance).
grid_data <- data.frame(
  x1 = runif(n),
  x2 = runif(n),
  Mean = runif(n),
  Upper = runif(n, 0.5, 1),
  Lower = runif(n, 0, 0.5),
  Variance = runif(n, 0, 0.1),
  class = max.col(Y_multi),
  max_prob = apply(Y_multi, 1, max)
)


test_that("my_2D_plot_fun generates plot without error", {
  # This test verifies that the function for binary/continuous plots runs correctly.
  expect_no_error({
    my_2D_plot_fun(
      var = "Mean",
      title = "Binary Predictive Mean Test",
      data = grid_data,
      X = X,
      y = y_bin
    )
  })
})


test_that("my_2D_plot_fun_class generates plot without error", {
  # This test verifies that the function for multiclass plots runs correctly.
  # It tests both classification and probability plotting modes.
  expect_no_error({
    # Test with multiclass classification data (classification = TRUE)
    my_2D_plot_fun_class(
      var = "class",
      title = "Multiclass Classification Test",
      data = grid_data,
      X = X,
      Y = Y_multi,
      classification = TRUE
    )
  })

  expect_no_error({
    # Test with multiclass probability data (classification = FALSE)
    my_2D_plot_fun_class(
      var = "max_prob",
      title = "Multiclass Max Probability Test",
      data = grid_data,
      X = X,
      Y = Y_multi,
      classification = FALSE
    )
  })
})

test_that("posterior_summary returns expected summary matrix", {
  ps <- posterior_summary(mean_vals = c(0.1, 0.2, 0.3), var_vals = c(0.01, 0.02, 0.03))
  expect_true(is.matrix(ps))
  expect_equal(rownames(ps), c("Posterior means", "Posterior variances"))
})

test_that(".make_plot_grid fixes non-displayed dimensions at training medians", {
  X <- matrix(c(1, 10, 100,
                2, 20, 200,
                3, 30, 300), ncol = 3, byrow = TRUE)
  grid <- expand.grid(x1 = c(0.1, 0.9), x3 = c(0.2, 0.8))

  Xnew <- .make_plot_grid(X, dims = c(1, 3), grid = grid)

  expect_equal(Xnew[, 2], rep(stats::median(X[, 2]), nrow(grid)))
  expect_equal(Xnew[, c(1, 3)], unname(as.matrix(grid)))
})

test_that("C++ Shepard interpolation matches reference R implementation", {
  Xtrain <- matrix(c(0, 0,
                     1, 0,
                     0, 1,
                     1, 1), ncol = 2, byrow = TRUE)
  Xquery <- matrix(c(0.5, 0.5,
                     0, 0), ncol = 2, byrow = TRUE)
  m <- c(10, 20, 30, 40)

  ref <- function(Xquery, Xtrain, m, power = 2) {
    vapply(seq_len(nrow(Xquery)), function(i) {
      dist_sq <- rowSums((t(t(Xtrain) - Xquery[i, ]))^2)
      if (any(dist_sq == 0)) {
        return(m[which(dist_sq == 0)[1]])
      }
      w <- dist_sq^(-power / 2)
      sum(w * m) / sum(w)
    }, numeric(1))
  }

  expect_equal(.bkp_shepard_m(Xquery, Xtrain, m), ref(Xquery, Xtrain, m))
})

test_that("C++ leave-one-out Shepard interpolation matches reference R implementation", {
  X <- matrix(c(0, 0,
                1, 0,
                0, 1,
                1, 1), ncol = 2, byrow = TRUE)
  m <- c(10, 20, 30, 40)

  ref <- vapply(seq_len(nrow(X)), function(i) {
    idx <- setdiff(seq_len(nrow(X)), i)
    dist_sq <- rowSums((t(t(X[idx, , drop = FALSE]) - X[i, ]))^2)
    w <- dist_sq^(-1)
    sum(w * m[idx]) / sum(w)
  }, numeric(1))

  expect_equal(.bkp_shepard_m_loo(X, m), ref)
})
