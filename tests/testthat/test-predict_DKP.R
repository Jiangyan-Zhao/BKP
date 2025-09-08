# This test file validates the functionality of the `predict` S3 method for
# DKP objects. It ensures that the function generates predictions without
# errors for both training and new data, and that the output structure
# and dimensions are correct.

test_that("predict.DKP generates predictions without errors for various priors", {
  # Set a seed for reproducibility to ensure consistent data
  set.seed(42)

  # -------------------------------------------------------------------------
  # Test Case 1: Noninformative Prior
  # -------------------------------------------------------------------------

  # Define true class probability function (3-class)
  true_pi_fun <- function(X) {
    p1 <- 1 / (1 + exp(-3 * X))
    p2 <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
    p_sum <- p1 + p2
    # Ensure probabilities sum to 1 and are valid
    p1_new <- p1 / (p_sum + 1)
    p2_new <- p2 / (p_sum + 1)
    p3_new <- 1 - p1_new - p2_new
    return(matrix(c(p1_new, p2_new, p3_new), ncol = 3))
  }

  d <- 1
  n <- 50
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- matrix(runif(n * d, -2, 2), n, d)
  true_pi <- true_pi_fun(X)
  m <- sample(150, n, replace = TRUE)
  Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
  colnames(Y) <- paste0("class", 1:3)

  # Fit DKP model with noninformative prior
  model_noninformative <- fit_DKP(X, Y, Xbounds = Xbounds, prior = "noninformative")

  # Define new prediction locations
  Xnew <- matrix(seq(-2, 2, length.out = 100), ncol = 1)

  # Predict on training data
  pred_train_noninformative <- predict(model_noninformative)
  expect_no_error(predict(model_noninformative))
  expect_identical(pred_train_noninformative$mean, fitted(model_noninformative))

  # Predict on new data
  pred_new_noninformative <- predict(model_noninformative, Xnew)
  expect_no_error(predict(model_noninformative, Xnew))
  expect_equal(nrow(pred_new_noninformative$mean), nrow(Xnew))
  expect_equal(ncol(pred_new_noninformative$mean), ncol(Y))
  expect_equal(nrow(pred_new_noninformative$variance), nrow(Xnew))
  expect_equal(ncol(pred_new_noninformative$variance), ncol(Y))
  expect_equal(nrow(pred_new_noninformative$lower), nrow(Xnew))
  expect_equal(ncol(pred_new_noninformative$lower), ncol(Y))
  expect_equal(nrow(pred_new_noninformative$upper), nrow(Xnew))
  expect_equal(ncol(pred_new_noninformative$upper), ncol(Y))

  # -------------------------------------------------------------------------
  # Test Case 2: Fixed Prior
  # -------------------------------------------------------------------------

  p0 <- c(0.2, 0.3, 0.5)
  model_fixed <- fit_DKP(X, Y, Xbounds = Xbounds, prior = "fixed", r0 = 10, p0 = p0)

  pred_train_fixed <- predict(model_fixed)
  expect_no_error(predict(model_fixed))
  expect_identical(pred_train_fixed$mean, fitted(model_fixed))

  pred_new_fixed <- predict(model_fixed, Xnew)
  expect_no_error(predict(model_fixed, Xnew))
  expect_equal(nrow(pred_new_fixed$mean), nrow(Xnew))
  expect_equal(ncol(pred_new_fixed$mean), ncol(Y))

  # -------------------------------------------------------------------------
  # Test Case 3: Adaptive Prior
  # -------------------------------------------------------------------------

  model_adaptive <- fit_DKP(X, Y, Xbounds = Xbounds, prior = "adaptive")

  pred_train_adaptive <- predict(model_adaptive)
  expect_no_error(predict(model_adaptive))
  expect_identical(pred_train_adaptive$mean, fitted(model_adaptive))

  pred_new_adaptive <- predict(model_adaptive, Xnew)
  expect_no_error(predict(model_adaptive, Xnew))
  expect_equal(nrow(pred_new_adaptive$mean), nrow(Xnew))
  expect_equal(ncol(pred_new_adaptive$mean), ncol(Y))
})
