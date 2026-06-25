test_that("plot.TwinBKP supports 1D base plots", {
  model <- make_twinbkp_model_1d()$model

  expect_silent(plot(model, n_grid = 10, engine = "base"))
  expect_silent(plot(model, n_grid = 10, engine = "base", show_global = FALSE))
})

test_that("plot.TwinBKP supports 2D base and ggplot plots and validates dims", {
  model_2d <- make_twinbkp_model_2d()$model

  expect_no_error(plot(model_2d, n_grid = 8, only_mean = TRUE, engine = "base"))
  expect_no_error(plot(model_2d, n_grid = 8, only_mean = TRUE, engine = "base", show_global = FALSE))
  expect_no_error(plot(model_2d, n_grid = 8, only_mean = TRUE, engine = "ggplot"))
  expect_no_error(plot(model_2d, n_grid = 8, only_mean = TRUE, engine = "ggplot", show_global = FALSE))
  expect_error(plot(model_2d, dims = c(1, 3), engine = "base"))
})
