test_that("ess = 'none' matches the standard BKP posterior update", {
  X <- matrix(c(0.05, 0.20, 0.45, 0.70, 0.90), ncol = 1)
  m <- c(10, 20, 15, 25, 30)
  y <- c(2, 12, 7, 15, 21)

  model <- fit_BKP(X, y, m, theta = 0.3, ess = "none")
  K <- kernel_matrix(model$Xnorm, theta = model$theta_opt,
                     kernel = model$kernel, isotropic = model$isotropic)
  prior_par <- get_prior(prior = model$prior, model = "BKP",
                         r0 = model$r0, p0 = model$p0,
                         y = model$y, m = model$m, K = K)

  expect_equal(model$ess, "none")
  expect_equal(model$alpha_n, prior_par$alpha0 + as.vector(K %*% model$y))
  expect_equal(model$beta_n, prior_par$beta0 + as.vector(K %*% (model$m - model$y)))
})

test_that("Shepard ESS targets observed trial sizes at training points", {
  X <- matrix(c(0.02, 0.18, 0.44, 0.73, 0.96), ncol = 1)
  m <- c(8, 17, 11, 23, 31)
  y <- c(3, 8, 4, 15, 22)

  model <- fit_BKP(X, y, m, theta = 0.25, ess = "shepard")
  data_size <- model$alpha_n + model$beta_n - model$alpha0 - model$beta0

  expect_equal(as.vector(data_size), as.numeric(m), tolerance = 1e-10)
  expect_equal(model$ess_info$m_shepard, as.numeric(m), tolerance = 1e-10)
})

test_that("Shepard ESS preserves the kernel-weighted empirical success proportion", {
  X <- matrix(c(0.10, 0.15, 0.35, 0.60, 0.85, 0.95), ncol = 2, byrow = TRUE)
  m <- c(10, 12, 18)
  y <- c(4, 9, 7)

  model <- fit_BKP(X, y, m, theta = 0.4, ess = "shepard")
  K <- kernel_matrix(model$Xnorm, theta = model$theta_opt,
                     kernel = model$kernel, isotropic = model$isotropic)
  raw_success <- as.vector(K %*% model$y)
  raw_trials <- as.vector(K %*% model$m)
  c_scale <- model$ess_info$scale

  positive <- raw_trials > 0
  expect_equal(
    raw_success[positive] / raw_trials[positive],
    (c_scale[positive] * raw_success[positive]) /
      (c_scale[positive] * raw_trials[positive]),
    tolerance = 1e-12
  )
})

test_that("predict works for Shepard ESS at new points without NA values", {
  X <- matrix(c(0.05, 0.25, 0.50, 0.80), ncol = 1)
  m <- c(10, 20, 15, 30)
  y <- c(2, 11, 8, 21)
  Xnew <- matrix(c(0.10, 0.40, 0.65, 0.90), ncol = 1)

  model <- fit_BKP(X, y, m, theta = 0.35, ess = "shepard")
  pred <- predict(model, Xnew = Xnew)

  expect_false(anyNA(pred$alpha_n))
  expect_false(anyNA(pred$beta_n))
  expect_false(anyNA(pred$mean))
  expect_true(all(is.finite(pred$alpha_n)))
  expect_true(all(is.finite(pred$beta_n)))
  expect_true(all(is.finite(pred$mean)))
})

test_that("Shepard ESS rejects duplicated input locations", {
  X <- matrix(c(0.10, 0.20, 0.10, 0.20, 0.75, 0.90), ncol = 2, byrow = TRUE)
  m <- c(10, 20, 15)
  y <- c(4, 9, 7)

  expect_error(
    fit_BKP(X, y, m, theta = 0.3, ess = "shepard"),
    "requires unique input locations"
  )
})
