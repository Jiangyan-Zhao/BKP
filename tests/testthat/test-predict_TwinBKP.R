test_that("predict.TwinBKP supports training, new-point, and count predictions", {
  set.seed(101)
  n <- 24
  X <- matrix(runif(n * 2), ncol = 2)
  m <- sample(8:20, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = 0.45)

  model <- fit_TwinBKP(X, y, m, theta_g = 0.3, theta_l = 0.4,
                       g = 8, runs = 2)

  pred <- predict(model)

  expect_s3_class(pred, "predict_TwinBKP")
  expect_s3_class(pred, "predict_BKP")
  expect_equal(pred$alpha_n, model$alpha_n)
  expect_equal(pred$beta_n, model$beta_n)
  expect_true(all(pred$mean >= 0 & pred$mean <= 1))

  Xnew <- matrix(runif(10), ncol = 2)
  pred <- predict(model, Xnew = Xnew)

  expect_equal(length(pred$mean), nrow(Xnew))
  expect_true(all(is.finite(pred$mean)))
  expect_true(all(is.finite(pred$lower)))
  expect_true(all(is.finite(pred$upper)))
  expect_true(all(pred$mean >= 0 & pred$mean <= 1))

  pred_count <- predict(model, Xnew = Xnew, type = "count",
                        Mnew = rep(20, nrow(Xnew)))

  expect_true("Mnew" %in% names(pred_count))
  expect_true(all(pred_count$mean >= 0))
  expect_true(all(pred_count$variance >= 0))
})

test_that("predict.TwinBKP supports classification summaries", {
  set.seed(102)
  n <- 20
  X <- matrix(runif(n * 2), ncol = 2)
  m <- rep(1, n)
  y <- rbinom(n, size = m, prob = 0.5)

  model <- fit_TwinBKP(X, y, m, theta_g = 0.3, theta_l = 0.4,
                       g = 6, runs = 2)

  pred <- predict(model)

  expect_true("class" %in% names(pred))
  expect_true(all(pred$class %in% c(0, 1)))
})

test_that("predict.TwinBKP validates inputs", {
  set.seed(103)
  n <- 20
  X <- matrix(runif(n * 2), ncol = 2)
  m <- sample(8:20, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = 0.45)
  Xnew <- matrix(runif(10), ncol = 2)

  model <- fit_TwinBKP(X, y, m, theta_g = 0.3, theta_l = 0.4,
                       g = 6, runs = 2)

  expect_error(predict(model, Xnew = matrix(runif(9), ncol = 3)),
               "The number of columns in 'Xnew' must match the original input dimension.")

  expect_error(predict(model, CI_level = 1),
               "'CI_level' must be a single numeric value strictly between 0 and 1.")

  expect_error(predict(model, threshold = 1),
               "'threshold' must be a single numeric value strictly between 0 and 1.")

  expect_error(predict(model, Xnew = Xnew, type = "count"),
               "When type = 'count' and Xnew is provided, 'Mnew' must also be provided.")
})

test_that("predict.TwinBKP supports Shepard ESS", {
  set.seed(104)
  n <- 20
  X <- matrix(runif(n * 2), ncol = 2)
  m <- rep(1, n)
  y <- rbinom(n, size = m, prob = 0.5)
  Xnew <- matrix(runif(10), ncol = 2)

  model <- fit_TwinBKP(X, y, m, theta_g = 0.3, theta_l = 0.4,
                       g = 6, runs = 2, ess = "shepard")

  pred_train <- predict(model)
  pred_new <- predict(model, Xnew = Xnew)

  expect_true(all(is.finite(pred_train$mean)))
  expect_true(all(is.finite(pred_new$mean)))
  expect_false(anyNA(pred_new$ess_info$scale))
})

test_that("predict.TwinBKP keeps adaptive prior positive when all kernel weights are zero", {
  set.seed(123)

  n <- 20
  X <- matrix(runif(n * 2), ncol = 2)
  m <- sample(8:20, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = 0.4)

  model <- fit_TwinBKP(
    X, y, m,
    prior = "adaptive",
    r0 = 2,
    global_kernel = "wendland",
    local_kernel = "wendland",
    theta_g = 0.05,
    theta_l = 0.05,
    g = 6,
    runs = 2
  )

  ## Far outside the normalized training domain; both Wendland components
  ## should have zero weight.
  Xnew <- matrix(c(10, 10), ncol = 2)

  pred <- predict(model, Xnew = Xnew)

  expect_true(all(is.finite(pred$alpha_n)))
  expect_true(all(is.finite(pred$beta_n)))
  expect_true(all(pred$alpha_n > 0))
  expect_true(all(pred$beta_n > 0))
  expect_true(all(is.finite(pred$mean)))
  expect_true(all(pred$mean >= 0 & pred$mean <= 1))
})
