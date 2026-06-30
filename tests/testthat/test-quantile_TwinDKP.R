test_that("quantile.TwinDKP returns DKP-style quantiles", {
  q <- quantile(make_twindkp_model_1d()$model, probs = c(0.025, 0.975))
  expect_equal(dim(q), c(40, 3, 2))
  }
)

test_that("quantile.TwinDKP handles single and invalid probabilities", {
  model <- make_twindkp_model_1d()$model

  q50 <- quantile(model, probs = 0.5)
  expect_equal(dim(q50), c(40, 3))

  expect_error(
    quantile(model, probs = numeric(0)),
    "nonempty finite numeric vector"
  )
})
