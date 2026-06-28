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


make_twindkp_model_2d_small <- function() {
  set.seed(2026)
  n <- 18
  X <- cbind(
    seq(0, 1, length.out = n),
    rep(c(0.2, 0.8), length.out = n)
  )

  eta1 <- 1.0 - X[, 1] + 0.5 * X[, 2]
  eta2 <- 0.2 + X[, 1] - 0.2 * X[, 2]
  eta3 <- -0.4 + 0.4 * X[, 1] + X[, 2]
  E <- cbind(eta1, eta2, eta3)
  P <- exp(E)
  P <- P / rowSums(P)

  m <- rep(8, n)
  Y <- t(vapply(
    seq_len(n),
    function(i) as.numeric(rmultinom(1, size = m[i], prob = P[i, ])),
    numeric(3)
  ))

  fit_TwinDKP(
    X, Y,
    prior = "fixed",
    r0 = 2,
    p0 = rep(1 / 3, 3),
    theta_g = 0.4,
    theta_l = 0.3,
    g = 6,
    l = 3,
    twins = 1,
    global_kernel = "gaussian",
    local_kernel = "wendland"
  )
}

test_that("plot.TwinDKP works for 2D base plots with global overlay", {
  fit <- make_twindkp_model_2d_small()
  expect_silent(plot(fit, n_grid = 5, engine = "base", show_global = TRUE))
})

test_that("plot.TwinDKP works for 2D ggplot plots with global overlay", {
  fit <- make_twindkp_model_2d_small()
  expect_silent(plot(fit, n_grid = 5, engine = "ggplot", show_global = TRUE))
})

test_that("plot.TwinDKP works for 2D single-label classification plots", {
  fit <- make_twindkp_model_2d_small()
  Y_onehot <- matrix(0, nrow = nrow(fit$Y), ncol = ncol(fit$Y))
  cls <- max.col(fit$Y)
  Y_onehot[cbind(seq_len(nrow(Y_onehot)), cls)] <- 1

  fit_cls <- fit_TwinDKP(
    fit$X, Y_onehot,
    prior = "fixed",
    r0 = 2,
    p0 = rep(1 / ncol(Y_onehot), ncol(Y_onehot)),
    theta_g = 0.4,
    theta_l = 0.3,
    g = 6,
    l = 3,
    twins = 1,
    global_kernel = "gaussian",
    local_kernel = "wendland"
  )

  expect_silent(plot(fit_cls, n_grid = 5, engine = "ggplot", show_global = TRUE))
})
