test_that("simulate.TwinDKP returns valid simulations", {
  s <- simulate(make_twindkp_model_1d()$model, nsim = 2, seed = 1)

  expect_s3_class(s, "simulate_TwinDKP")
  expect_equal(dim(s$samples), c(40, 3, 2))
  expect_null(s$class)
})

test_that("simulate.TwinDKP returns MAP classes for single-label data", {
  fx <- make_twindkp_data()
  Y_onehot <- matrix(0, nrow = nrow(fx$Y), ncol = ncol(fx$Y))
  cls <- max.col(fx$Y)
  Y_onehot[cbind(seq_len(nrow(Y_onehot)), cls)] <- 1

  fit <- fit_TwinDKP(
    fx$X, Y_onehot,
    prior = "fixed",
    r0 = 2,
    p0 = rep(1 / ncol(Y_onehot), ncol(Y_onehot)),
    theta_g = 0.4,
    theta_l = 0.3,
    g = 10,
    l = 5,
    twins = 1,
    global_kernel = "gaussian",
    local_kernel = "wendland"
  )

  s <- simulate(fit, nsim = 2, seed = 1)

  expect_s3_class(s, "simulate_TwinDKP")
  expect_equal(dim(s$samples), c(nrow(fx$X), ncol(Y_onehot), 2))
  expect_equal(dim(s$class), c(nrow(fx$X), 2))
  expect_true(all(s$class %in% seq_len(ncol(Y_onehot))))
})
