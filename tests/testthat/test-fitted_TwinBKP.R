test_that("fitted.TwinBKP returns training posterior means", {
  fit <- make_twinbkp_model_1d()
  model <- fit$model
  fv <- fitted(model)

  expect_type(fv, "double")
  expect_length(fv, nrow(fit$X))
  expect_true(all(fv >= 0 & fv <= 1))
  expect_equal(fv, model$alpha_n / (model$alpha_n + model$beta_n))
  expect_equal(fv, predict(model)$mean)
})
