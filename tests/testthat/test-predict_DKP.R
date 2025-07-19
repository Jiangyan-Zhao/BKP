# -------------------------- Test Setup ---------------------------
set.seed(123)
n_model <- 30
d_model <- 2
q_model <- 3 # Number of classes
X_model <- matrix(runif(n_model * d_model, 0, 1), n_model, d_model)
m_sum_model <- sample(50, n_model, replace = TRUE)
Y_model <- t(sapply(1:n_model, function(i) rmultinom(1, size = m_sum_model[i], prob = rep(1/q_model, q_model))))

# Create a dummy DKP model object (minimal for testing predict.DKP)
dummy_DKP_model <- list(
  theta_opt = rep(1, d_model),
  kernel = "gaussian",
  prior = "noninformative",
  r0 = 2,
  p0 = rep(1/q_model, q_model),
  X = X_model,
  Xnorm = X_model, # Assuming X is already normalized for simplicity
  Xbounds = cbind(rep(0, d_model), rep(1, d_model)),
  Y = Y_model
)
class(dummy_DKP_model) <- "DKP"

n_new_test <- 5
Xnew_test <- matrix(runif(n_new_test * d_model, 0, 1), n_new_test, d_model)

# -------------------------- Test Context: Basic Functionality ---------------------------
test_that("Basic Functionality: predict.DKP returns a list with expected components", {
  prediction <- predict.DKP(dummy_DKP_model, Xnew_test)
  expect_type(prediction, "list")
  expect_named(prediction, c("Xnew", "mean", "variance", "lower", "upper", "CI_level"))
})

test_that("predicted mean, variance, lower, upper have correct dimensions", {
  prediction <- predict.DKP(dummy_DKP_model, Xnew_test)
  expect_equal(dim(prediction$mean), c(n_new_test, q_model))
  expect_equal(dim(prediction$variance), c(n_new_test, q_model))
  expect_equal(dim(prediction$lower), c(n_new_test, q_model))
  expect_equal(dim(prediction$upper), c(n_new_test, q_model))
})

test_that("CI_level is correctly passed through", {
  prediction <- predict.DKP(dummy_DKP_model, Xnew_test, CI_level = 0.1)
  expect_equal(prediction$CI_level, 0.1)
})

test_that("'class' output is correct for classification data (rowSums(Y)==1)", {
  # Create a dummy DKP model for classification data (each row in Y sums to 1)
  dummy_DKP_model_class <- dummy_DKP_model
  dummy_DKP_model_class$Y <- t(sapply(1:n_model, function(i) {
    p <- rep(0, q_model)
    p[sample(1:q_model, 1)] <- 1 # Only one class is 1, others 0
    p
  }))

  # Mock kernel_matrix and get_prior_dkp to control alpha_n, thus pi_mean
  mock_kernel_matrix <- function(Xnew_norm, Xnorm, theta, kernel) matrix(1, nrow(Xnew_norm), nrow(Xnorm))
  mock_get_prior_dkp <- function(prior, r0, p0, Y, K) matrix(c(10,1,1, 1,10,1, 1,1,10, 10,1,1, 1,10,1), byrow=TRUE, nrow=n_new_test, ncol=q_model)

  with_mocked_bindings(
    {
      prediction <- predict.DKP(dummy_DKP_model_class, Xnew_test)
      # Based on mock_get_prior_dkp, alpha_n will be (10+1, 1+1, 1+1), (1+1, 10+1, 1+1)...
      # Expected classes will be 1, 2, 3, 1, 2 for the 5 rows
      expect_equal(as.vector(prediction$class), c(1, 2, 3, 1, 2))
    },
    kernel_matrix = mock_kernel_matrix,
    get_prior_dkp = mock_get_prior_dkp
  )
})

# -------------------------- Test Context: Input Validation ---------------------------

test_that("Input Validation: input 'object' must be of class 'DKP'", {
  expect_error(predict.DKP(list(), Xnew_test))
  expect_error(predict.DKP(matrix(), Xnew_test))
})

test_that("Input Validation: 'Xnew' must have the same number of columns as original X", {
  Xnew_wrong_dim <- matrix(runif(n_new_test * (d_model + 1)), n_new_test, d_model + 1)
  expect_error(predict.DKP(dummy_DKP_model, Xnew_wrong_dim))
})

test_that("Input Validation: vector 'Xnew' input is handled correctly", {
  Xnew_vector <- Xnew_test[1, ] # Take first row as a vector
  prediction <- predict.DKP(dummy_DKP_model, Xnew_vector)
  expect_equal(dim(prediction$mean), c(1, q_model))
  expect_equal(prediction$Xnew, matrix(Xnew_vector, nrow = 1))
})
