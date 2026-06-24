test_that("quantile.TwinBKP returns posterior quantiles and validates probs", {
  fit <- make_twinbkp_model_1d()
  model <- fit$model

  q <- quantile(model)
  expect_true(is.matrix(q))
  expect_equal(dim(q), c(nrow(fit$X), 3L))
  q50 <- quantile(model, probs = 0.5)
  expect_type(q50, "double")
  expect_length(q50, nrow(fit$X))
  expect_true(all(q >= 0 & q <= 1))
  expect_true(all(q50 >= 0 & q50 <= 1))
  expect_error(quantile(model, probs = -0.1))
  expect_error(quantile(model, probs = 1.1))
})
