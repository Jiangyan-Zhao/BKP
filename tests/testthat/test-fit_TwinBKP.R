test_that("TwinBKP ess = 'none' preserves default fit and prediction behavior", {
  set.seed(301)
  X <- matrix(seq(0.05, 0.95, length.out = 12), ncol = 1)
  m <- c(10, 12, 9, 11, 13, 8, 10, 12, 9, 11, 13, 8)
  y <- c(2, 3, 4, 5, 6, 3, 4, 7, 5, 6, 8, 4)

  default_model <- fit_TwinBKP(X, y, m, theta = 0.4, g_nums = 4)
  none_model <- fit_TwinBKP(X, y, m, theta = 0.4, g_nums = 4, ess = "none")

  expect_equal(none_model$ess, "none")
  expect_equal(none_model$ess_info_global$scale, rep(1, length(none_model$global_idx)))
  expect_equal(none_model$alpha_n_global, default_model$alpha_n_global)
  expect_equal(none_model$beta_n_global, default_model$beta_n_global)

  Xnew <- matrix(c(0.2, 0.55, 0.9), ncol = 1)
  set.seed(302)
  default_pred <- predict(default_model, Xnew = Xnew, l_nums = 4, v_nums = 4)
  set.seed(302)
  none_pred <- predict(none_model, Xnew = Xnew, l_nums = 4, v_nums = 4)

  expect_equal(none_pred$alpha_n, default_pred$alpha_n)
  expect_equal(none_pred$beta_n, default_pred$beta_n)
  expect_equal(none_pred$mean, default_pred$mean)
  expect_equal(none_pred$variance, default_pred$variance)
})

test_that("TwinBKP supports Shepard ESS at fit and prediction time", {
  set.seed(303)
  X <- matrix(seq(0.04, 0.96, length.out = 14), ncol = 1)
  m <- c(8, 10, 12, 9, 11, 13, 10, 8, 12, 14, 9, 11, 13, 10)
  y <- c(1, 3, 5, 4, 6, 7, 4, 3, 6, 9, 5, 7, 8, 6)

  model <- fit_TwinBKP(X, y, m, theta = 0.35, g_nums = 5, ess = "shepard")

  expect_equal(model$ess, "shepard")
  expect_named(model$ess_info_global, c("scale", "m_kernel", "m_shepard", "m_target", "rho"))
  expect_false(anyNA(model$alpha_n_global))
  expect_false(anyNA(model$beta_n_global))

  Xnew <- rbind(X[model$global_idx[1], , drop = FALSE], matrix(c(0.18, 0.73), ncol = 1))
  set.seed(304)
  pred <- predict(model, Xnew = Xnew, l_nums = 5, v_nums = 5)

  expect_equal(pred$ess, "shepard")
  expect_false(anyNA(pred$mean))
  expect_false(anyNA(pred$variance))
  expect_false(anyNA(pred$lower))
  expect_false(anyNA(pred$upper))

  exact_ess <- pred$ess_info[[1]]
  calibrated_ess <- exact_ess$scale * exact_ess$m_kernel
  expect_equal(as.numeric(calibrated_ess), as.numeric(m[model$global_idx[1]]), tolerance = 1e-8)
})
