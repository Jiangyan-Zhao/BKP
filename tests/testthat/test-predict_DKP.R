test_that("predict.DKP returns expected structure and values", {
  # 0. Define true class probability function (3-class)
  true_pi_fun <- function(X) {
    p <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
    return(matrix(c(p/2, p/2, 1 - p), nrow = length(p)))
  }

  # 1. Simulate training data
  n <- 30
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- tgp::lhs(n = n, rect = Xbounds)
  true_pi <- true_pi_fun(X)
  m <- sample(100, n, replace = TRUE)
  Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))

  # 2. Fit DKP model
  model <- fit.DKP(X, Y, Xbounds = Xbounds)

  # 3. Predict on new input
  n_Xnew <- 10
  Xnew <- matrix(seq(-2, 2, length.out = n_Xnew), ncol = 1)
  prediction <- predict(model, Xnew)

  # 4. Check structure
  expect_type(prediction, "list")
  expect_in(names(prediction),
            c("Xnew", "mean", "variance", "lower", "upper", "CI_level", "class"))

  # 5. Dimensions and types
  expect_equal(nrow(prediction$Xnew), n_Xnew)
  expect_equal(dim(prediction$mean), c(n_Xnew, 3))
  expect_equal(dim(prediction$variance), c(n_Xnew, 3))
  expect_equal(dim(prediction$lower), c(n_Xnew, 3))
  expect_equal(dim(prediction$upper), c(n_Xnew, 3))
  expect_true(all(prediction$class %in% 1:3))

  # 6. Prediction values should be probabilities
  expect_true(all(prediction$mean >= 0 & prediction$mean <= 1))
  expect_true(all(prediction$lower >= 0 & prediction$upper <= 1))
  expect_true(all(prediction$upper >= 0 & prediction$upper <= 1))
})
