# tests/testthat/test-loss_fun_dkp.R

library(testthat)


# -------------------------- Test Setup ---------------------------
set.seed(42)
n_test <- 10
d_test <- 2
q_test <- 3 # Number of classes
Xnorm_test <- matrix(runif(n_test * d_test), n_test, d_test)
# Generate multinomial counts for Y
m_sum_test <- sample(50, n_test, replace = TRUE)
Y_test <- t(sapply(1:n_test, function(i) rmultinom(1, size = m_sum_test[i], prob = rep(1/q_test, q_test))))
gamma_test <- rep(0.5, d_test) # log10(theta) = 0.5 -> theta = sqrt(10)

# -------------------------- Test Context: Basic Functionality ---------------------------

test_that("Basic Functionality: loss_fun_dkp returns a single double value", {
  loss_val <- loss_fun_dkp(gamma_test, Xnorm_test, Y_test)
  expect_type(loss_val, "double")
  expect_length(loss_val, 1)
})

test_that("different prior types are accepted", {
  expect_silent(loss_fun_dkp(gamma_test, Xnorm_test, Y_test, prior = "noninformative"))
  expect_silent(loss_fun_dkp(gamma_test, Xnorm_test, Y_test, prior = "fixed", r0 = 5, p0 = rep(1/3, q_test)))
  expect_silent(loss_fun_dkp(gamma_test, Xnorm_test, Y_test, prior = "adaptive", r0 = 10))
})

test_that("different kernel types are accepted", {
  expect_silent(loss_fun_dkp(gamma_test, Xnorm_test, Y_test, kernel = "gaussian"))
  expect_silent(loss_fun_dkp(gamma_test, Xnorm_test, Y_test, kernel = "matern52"))
  expect_silent(loss_fun_dkp(gamma_test, Xnorm_test, Y_test, kernel = "matern32"))
})

# -------------------------- Test Context: Input Validation ---------------------------

test_that("Input Validation: unsupported loss type throws error", {
  expect_error(loss_fun_dkp(gamma_test, Xnorm_test, Y_test, loss = "invalid_loss"))
})

test_that("Input Validation: unsupported prior type throws error", {
  expect_error(loss_fun_dkp(gamma_test, Xnorm_test, Y_test, prior = "invalid_prior"))
})

test_that("Input Validation: unsupported kernel type throws error", {
  expect_error(loss_fun_dkp(gamma_test, Xnorm_test, Y_test, kernel = "invalid_kernel"))
})

test_that("Input Validation: numerical stabilization prevents NaNs/Infs in pi_hat", {
  # Create a scenario where alpha_n could be very small
  mock_kernel_matrix_zero <- function(Xnorm, theta, kernel) matrix(0, nrow(Xnorm), nrow(Xnorm))
  mock_get_prior_small <- function(prior, r0, p0, Y, K) matrix(1e-15, nrow(Y), ncol(Y))

  with_mocked_bindings(
    {
      loss_val <- loss_fun_dkp(gamma_test, Xnorm_test, Y_test)
      expect_false(is.nan(loss_val))
      expect_false(is.infinite(loss_val))
    },
    kernel_matrix = mock_kernel_matrix_zero,
    get_prior_dkp = mock_get_prior_small
  )
})
