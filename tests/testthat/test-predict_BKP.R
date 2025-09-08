# This test file validates the functionality of the `predict` S3 method for
# BKP objects. It ensures that the function generates predictions without
# errors for both training and new data, and that the output structure
# and dimensions are correct.

test_that("predict.BKP generates predictions without errors for various priors", {
  # Set a seed for reproducibility to ensure consistent data
  set.seed(42)

  # -------------------------------------------------------------------------
  # Test Case 1: Noninformative Prior
  # -------------------------------------------------------------------------

  # Define true success probability function
  true_pi_fun <- function(x) {
    (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
  }

  d <- 1
  n <- 50
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- matrix(runif(n * d, -2, 2), n, d)
  true_pi <- true_pi_fun(X)
  m <- sample(100, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = true_pi)

  # Fit BKP model with noninformative prior
  model_noninformative <- fit_BKP(X, y, m, prior = "noninformative", Xbounds = Xbounds)

  # Define new prediction locations
  Xnew <- matrix(seq(-2, 2, length.out = 100), ncol = 1)

  # Predict on training data
  pred_train_noninformative <- predict(model_noninformative)
  expect_no_error(predict(model_noninformative))
  expect_identical(pred_train_noninformative$mean, fitted(model_noninformative))

  # Predict on new data
  pred_new_noninformative <- predict(model_noninformative, Xnew)
  expect_no_error(predict(model_noninformative, Xnew))
  expect_equal(length(pred_new_noninformative$mean), nrow(Xnew))
  expect_equal(length(pred_new_noninformative$variance), nrow(Xnew))
  expect_equal(length(pred_new_noninformative$lower), nrow(Xnew))
  expect_equal(length(pred_new_noninformative$upper), nrow(Xnew))
  expect_true(all(pred_new_noninformative$lower <= pred_new_noninformative$upper))

  # -------------------------------------------------------------------------
  # Test Case 2: Fixed Prior
  # -------------------------------------------------------------------------

  model_fixed <- fit_BKP(X, y, m, prior = "fixed", r0 = 10, p0 = 0.5, Xbounds = Xbounds)

  pred_train_fixed <- predict(model_fixed)
  expect_no_error(predict(model_fixed))
  expect_identical(pred_train_fixed$mean, fitted(model_fixed))

  pred_new_fixed <- predict(model_fixed, Xnew)
  expect_no_error(predict(model_fixed, Xnew))
  expect_equal(length(pred_new_fixed$mean), nrow(Xnew))
  expect_true(all(pred_new_fixed$lower <= pred_new_fixed$upper))

  # -------------------------------------------------------------------------
  # Test Case 3: Adaptive Prior
  # -------------------------------------------------------------------------

  model_adaptive <- fit_BKP(X, y, m, prior = "adaptive", Xbounds = Xbounds)

  pred_train_adaptive <- predict(model_adaptive)
  expect_no_error(predict(model_adaptive))
  expect_identical(pred_train_adaptive$mean, fitted(model_adaptive))

  pred_new_adaptive <- predict(model_adaptive, Xnew)
  expect_no_error(predict(model_adaptive, Xnew))
  expect_equal(length(pred_new_adaptive$mean), nrow(Xnew))
  expect_true(all(pred_new_adaptive$lower <= pred_new_adaptive$upper))
})
