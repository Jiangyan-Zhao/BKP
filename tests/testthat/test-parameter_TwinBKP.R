test_that("parameter.TwinBKP returns key model parameters", {
  model <- make_twinbkp_model_1d()$model
  pars <- parameter(model)

  expect_type(pars, "list")
  expect_true(all(c(
    "theta", "theta_g", "theta_l", "global_kernel", "local_kernel",
    "alpha_n", "beta_n", "global_indices", "control"
  ) %in% names(pars)))
  expect_equal(pars$theta, model$theta_opt)
  expect_equal(pars$theta_g, model$theta_g)
  expect_equal(pars$theta_l, model$theta_l)
  expect_equal(pars$alpha_n, model$alpha_n)
  expect_equal(pars$beta_n, model$beta_n)
})
