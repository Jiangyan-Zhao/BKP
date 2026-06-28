test_that("plot.TwinDKP works for 1D base plots", {
  fit <- make_twindkp_model_1d()$model
  expect_silent(plot(fit, n_grid = 5, engine = "base"))
})

test_that("plot.TwinDKP works for 1D ggplot plots", {
  fit <- make_twindkp_model_1d()$model
  expect_silent(plot(fit, n_grid = 5, engine = "ggplot"))
})

test_that("plot.TwinDKP validates plotting arguments", {
  fit <- make_twindkp_model_1d()$model
  expect_error(plot(fit, n_grid = 0), "n_grid")
  expect_error(plot(fit, only_mean = c(TRUE, FALSE)), "only_mean")
  expect_error(plot(fit, show_global = c(TRUE, FALSE)), "show_global")
})
