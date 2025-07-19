# -------------------------- Test Setup ---------------------------
set.seed(42)
n_test <- 10
d_test <- 2
Xnorm_test <- matrix(runif(n_test * d_test), n_test, d_test)
m_test <- rep(10, n_test)
y_test <- rbinom(n_test, size = m_test, prob = runif(n_test, 0.1, 0.9))
gamma_test <- rep(0.5, d_test) # log10(theta) = 0.5 -> theta = sqrt(10)

# -------------------------- Test Context: Basic Functionality ---------------------------

test_that("Basic Functionality: loss_fun returns a single double value", {
  loss_val <- loss_fun(gamma_test, Xnorm_test, y_test, m_test)
  expect_type(loss_val, "double")
  expect_length(loss_val, 1)
})


test_that("different prior types are accepted", {
  expect_silent(loss_fun(gamma_test, Xnorm_test, y_test, m_test, prior = "noninformative"))
  expect_silent(loss_fun(gamma_test, Xnorm_test, y_test, m_test, prior = "fixed", r0 = 5, p0 = 0.5))
  expect_silent(loss_fun(gamma_test, Xnorm_test, y_test, m_test, prior = "adaptive", r0 = 10))
})

test_that("different kernel types are accepted", {
  expect_silent(loss_fun(gamma_test, Xnorm_test, y_test, m_test, kernel = "gaussian"))
  expect_silent(loss_fun(gamma_test, Xnorm_test, y_test, m_test, kernel = "matern52"))
  expect_silent(loss_fun(gamma_test, Xnorm_test, y_test, m_test, kernel = "matern32"))
})

# -------------------------- Test Context: Input Validation ---------------------------

test_that("Input Validation: unsupported loss type throws error", {
  expect_error(loss_fun(gamma_test, Xnorm_test, y_test, m_test, loss = "invalid_loss"))
})

test_that("Input Validation: unsupported prior type throws error", {
  expect_error(loss_fun(gamma_test, Xnorm_test, y_test, m_test, prior = "invalid_prior"))
})

test_that("Input Validation: unsupported kernel type throws error", {
  expect_error(loss_fun(gamma_test, Xnorm_test, y_test, m_test, kernel = "invalid_kernel"))
})

test_that("Input Validation: numerical stabilization prevents NaNs/Infs in pi_hat", {
  # Create a scenario where alpha_n or beta_n could be very small
  # Mock prior and kernel_matrix to force small alpha/beta
  mock_kernel_matrix_zero <- function(Xnorm, theta, kernel) matrix(0, nrow(Xnorm), nrow(Xnorm))
  mock_get_prior_small <- function(prior, r0, p0, y, m, K) list(alpha0 = rep(1e-15, length(y)), beta0 = rep(1e-15, length(y)))

  with_mocked_bindings(
    {
      loss_val <- loss_fun(gamma_test, Xnorm_test, y_test, m_test)
      expect_false(is.nan(loss_val))
      expect_false(is.infinite(loss_val))
    },
    kernel_matrix = mock_kernel_matrix_zero,
    get_prior = mock_get_prior_small
  )
})
