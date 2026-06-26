test_that("fit_TwinDKP returns expected structure", {
  fx <- make_twindkp_model_1d(); fit <- fx$model
  expect_s3_class(fit, "TwinDKP")
  expect_length(fit$global_indices, fit$control$g)
  expect_equal(fit$control$g_target, 10)
  expect_equal(dim(fit$prob), c(nrow(fx$X), ncol(fx$Y)))
  expect_equal(as.numeric(rowSums(fit$prob)), rep(1, nrow(fx$X)), tolerance = 1e-8)
})
test_that("fit_TwinDKP validates inputs", {
  fx <- make_twindkp_data()
  expect_error(fit_TwinDKP(fx$X, fx$Y, prior = "fixed", p0 = c(.5, .5), theta_g = .4, theta_l = .3, g = 10, l = 5), "p0")
  badY <- fx$Y; badY[1, 1] <- -1
  expect_error(fit_TwinDKP(fx$X, badY, theta_g = .4, theta_l = .3, g = 10, l = 5), "nonnegative")
  expect_error(fit_TwinDKP(fx$X, fx$Y, theta_g = .4, theta_l = .3, g = 10, l = 5, control = list(ess = "shepard")), "does not support ESS")
})
