test_that("summary.TwinBKP returns expected fields and posterior summaries", {
  fit <- make_twinbkp_model_1d()
  model <- fit$model
  s <- summary(model)

  expect_type(s, "list")
  expect_s3_class(s, "summary_TwinBKP")
  expect_true(all(c(
    "n_obs", "input_dim", "kernel", "global_kernel", "local_kernel",
    "isotropic", "theta_opt", "theta_g", "theta_l", "loss", "loss_min",
    "prior", "r0", "p0", "global_size", "global_target", "local_size",
    "twins", "post_mean", "post_var"
  ) %in% names(s)))

  pred <- predict(model)
  expect_equal(s$post_mean, pred$mean)
  expect_equal(s$post_var, pred$variance)
  expect_output(print(s), "Twin Beta Kernel Process")
})
