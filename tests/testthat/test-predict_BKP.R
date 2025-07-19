test_that("predict.BKP returns correct structure and valid values", {
  # 0. Define true success probability function
  true_pi_fun <- function(x) {
    (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
  }

  # 1. Simulate binary outcome data
  set.seed(123)
  n <- 30
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- tgp::lhs(n = n, rect = Xbounds)
  true_pi <- true_pi_fun(X)
  m <- sample(100, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = true_pi)

  # 2. Fit BKP model
  model <- fit.BKP(X, y, m, Xbounds = Xbounds)

  # 3. Make predictions
  n_Xnew <- 10
  Xnew <- matrix(seq(-2, 2, length.out = n_Xnew), ncol = 1)
  prediction <- predict(model, Xnew)

  # 4. Check structure
  expect_type(prediction, "list")
  expect_in(names(prediction),
            c("Xnew", "mean", "variance", "lower", "upper", "CI_level", "class"))

  # 5. Check dimensions
  expect_equal(nrow(prediction$Xnew), n_Xnew)
  expect_length(prediction$mean, n_Xnew)
  expect_length(prediction$variance, n_Xnew)
  expect_length(prediction$lower, n_Xnew)
  expect_length(prediction$upper, n_Xnew)
  expect_true(all(prediction$class %in% c(0L, 1L)))

  # 5. Check probability bounds
  expect_true(all(prediction$mean >= 0 & prediction$mean <= 1))
  expect_true(all(prediction$lower >= 0 & prediction$lower <= 1))
  expect_true(all(prediction$upper >= 0 & prediction$upper <= 1))
})
