test_that("parameter.TwinDKP returns parameters", {
  p <- parameter(make_twindkp_model_1d()$model)
  expect_named(p, c(
    "theta", "theta_g", "theta_l", "alpha_n",
    "global_kernel", "local_kernel", "global_indices", "control"))
  }
)
