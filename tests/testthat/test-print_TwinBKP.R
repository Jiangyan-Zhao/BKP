test_that("TwinBKP print methods are silent", {
  model <- make_twinbkp_model_1d()$model

  expect_silent(print(model))
  expect_silent(print(predict(model)))
  expect_silent(print(simulate(model, nsim = 2, seed = 1)))
})
