# tests/testthat/test-predict_BKP.R

library(testthat)
# -------------------------- Test Setup ---------------------------
set.seed(123)
n_model <- 30
d_model <- 2
X_model <- matrix(runif(n_model * d_model, 0, 1), n_model, d_model)
m_model <- rep(10, n_model)
y_model <- rbinom(n_model, size = m_model, prob = runif(n_model, 0.1, 0.9))

# Create a dummy BKP model object (minimal for testing predict.BKP)
dummy_BKP_model <- list(
  theta_opt = rep(1, d_model),
  kernel = "gaussian",
  prior = "noninformative",
  r0 = 2,
  p0 = 0.5,
  X = X_model,
  Xnorm = X_model, # Assuming X is already normalized for simplicity
  Xbounds = cbind(rep(0, d_model), rep(1, d_model)),
  y = y_model,
  m = m_model
)
class(dummy_BKP_model) <- "BKP"

n_new_test <- 5
Xnew_test <- matrix(runif(n_new_test * d_model, 0, 1), n_new_test, d_model)

# -------------------------- Test Context: Basic Functionality ---------------------------
test_that("Basic Functionality: predict.BKP returns a list with expected components", {
  prediction <- predict.BKP(dummy_BKP_model, Xnew_test)
  expect_type(prediction, "list")
  expect_named(prediction, c("Xnew", "mean", "variance", "lower", "upper", "CI_level"))
})

test_that("predicted mean, variance, lower, upper have correct dimensions", {
  prediction <- predict.BKP(dummy_BKP_model, Xnew_test)
  expect_equal(length(prediction$mean), n_new_test)
  expect_equal(length(prediction$variance), n_new_test)
  expect_equal(length(prediction$lower), n_new_test)
  expect_equal(length(prediction$upper), n_new_test)
})

test_that("CI_level is correctly passed through", {
  prediction <- predict.BKP(dummy_BKP_model, Xnew_test, CI_level = 0.1)
  expect_equal(prediction$CI_level, 0.1)
})


# -------------------------- Test Context: Input Validation ---------------------------

test_that("Input Validation: input 'object' must be of class 'BKP'", {
  expect_error(predict.BKP(list(), Xnew_test)) # Not a BKP object
  expect_error(predict.BKP(matrix(), Xnew_test)) # Not a BKP object
})

test_that("Input Validation: 'Xnew' must have the same number of columns as original X", {
  Xnew_wrong_dim <- matrix(runif(n_new_test * (d_model + 1)), n_new_test, d_model + 1)
  expect_error(predict.BKP(dummy_BKP_model, Xnew_wrong_dim))
})

test_that("Input Validation: vector 'Xnew' input is handled correctly", {
  Xnew_vector <- Xnew_test[1, ] # Take first row as a vector
  prediction <- predict.BKP(dummy_BKP_model, Xnew_vector)
  expect_equal(length(prediction$mean), 1)
  expect_equal(prediction$Xnew, matrix(Xnew_vector, nrow = 1))
})
