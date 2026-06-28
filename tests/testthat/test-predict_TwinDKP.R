test_that("predict.TwinDKP returns probability and count predictions", {
  fx <- make_twindkp_model_1d()
  fit <- fx$model
  Xnew <- matrix(c(0.2, 0.5, 0.8), ncol = 1)

  p <- predict(fit, Xnew, type = "probability")
  expect_s3_class(p, "predict_TwinDKP")
  expect_equal(dim(p$mean), c(3, 3))
  expect_equal(as.numeric(rowSums(p$mean)), rep(1, 3), tolerance = 1e-8)
  expect_null(p$class)

  cnt <- predict(fit, Xnew, type = "count", Mnew = 12)
  expect_s3_class(cnt, "predict_TwinDKP")
  expect_equal(dim(cnt$mean), c(3, 3))
  expect_equal(dim(cnt$variance), c(3, 3))
  expect_equal(dim(cnt$lower), c(3, 3))
  expect_equal(dim(cnt$upper), c(3, 3))
  expect_equal(cnt$Mnew, rep(12L, 3))
  expect_null(cnt$class)
})

test_that("predict.TwinDKP validates Xnew", {
  fit <- make_twindkp_model_1d()$model
  expect_error(predict(fit, matrix(NA_real_, ncol = 1)), "finite")
  expect_error(predict(fit, matrix(NaN, ncol = 1)), "finite")
  expect_error(predict(fit, matrix(Inf, ncol = 1)), "finite")
})
