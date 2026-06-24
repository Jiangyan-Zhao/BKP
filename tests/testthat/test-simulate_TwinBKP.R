test_that("simulate.TwinBKP returns posterior simulations", {
  fit <- make_twinbkp_model_1d()
  model <- fit$model

  sim <- simulate(model, nsim = 3, seed = 1)
  expect_s3_class(sim, "simulate_TwinBKP")
  expect_equal(dim(sim$samples), c(nrow(fit$X), 3L))
  expect_length(sim$mean, nrow(fit$X))
  expect_true(all(sim$samples >= 0 & sim$samples <= 1))

  Xnew <- matrix(c(0.1, 0.5, 0.9), ncol = 1)
  sim_new <- simulate(model, Xnew = Xnew, nsim = 2, seed = 1)
  expect_equal(dim(sim_new$samples), c(nrow(Xnew), 2L))

  sim_class <- simulate(model, threshold = 0.5, nsim = 2, seed = 1)
  expect_true(all(sim_class$class %in% c(0L, 1L)))
})

test_that("simulate.TwinBKP validates arguments", {
  model <- make_twinbkp_model_1d()$model

  expect_error(simulate(model, nsim = 0))
  expect_error(simulate(model, seed = c(1, 2)))
  expect_error(simulate(model, Xnew = matrix(1, nrow = 1, ncol = 2)))
  expect_error(simulate(model, threshold = 1.1))
})
