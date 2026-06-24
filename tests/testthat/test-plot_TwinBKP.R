test_that("plot.TwinBKP supports 1D base plots", {
  model <- make_twinbkp_model_1d()$model

  expect_silent(plot(model, n_grid = 10, engine = "base"))
  expect_silent(plot(model, n_grid = 10, engine = "base", show_global = FALSE))
})

test_that("plot.TwinBKP supports 2D base plots and validates dims", {
  model <- make_twinbkp_model_2d()$model

  expect_silent(plot(model, n_grid = 8, only_mean = TRUE, engine = "base"))
  expect_error(plot(model, dims = c(1, 3), engine = "base"))
})
