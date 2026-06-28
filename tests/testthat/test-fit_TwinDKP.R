test_that("fit_TwinDKP returns expected structure", {
  fx <- make_twindkp_model_1d()
  fit <- fx$model
  expect_s3_class(fit, "TwinDKP")
  expect_length(fit$global_indices, fit$control$g)
  expect_equal(fit$control$g_target, 10)
  expect_equal(dim(fit$prob), c(nrow(fx$X), ncol(fx$Y)))
  expect_equal(as.numeric(rowSums(fit$prob)), rep(1, nrow(fx$X)), tolerance = 1e-8)
})
test_that("fit_TwinDKP validates inputs", {
  fx <- make_twindkp_data()
  expect_error(fit_TwinDKP(fx$X, fx$Y, prior = "fixed", p0 = c(.5, .5), theta_g = .4, theta_l = .3, g = 10, l = 5), "p0")
  badY <- fx$Y
  badY[1, 1] <- -1
  expect_error(fit_TwinDKP(fx$X, badY, theta_g = .4, theta_l = .3, g = 10, l = 5), "nonnegative")
  expect_error(
    fit_TwinDKP(
      fx$X, fx$Y, theta_g = .4, theta_l = .3, g = 10, l = 5,
      n_multi_start = 0
    ),
    "n_multi_start"
  )
  expect_error(
    fit_TwinDKP(
      fx$X, fx$Y, theta_g = .4, theta_l = .3, g = 10, l = 5,
      n_threads = 0
    ),
    "n_threads"
  )
})

.twindkp_dense_reference <- function(Xquery_norm, Xtrain_norm, Y, g_indices, local_indices,
                                     theta_g, theta_l, global_kernel, local_kernel,
                                     isotropic, prior, r0, p0) {
  t <- nrow(Xquery_norm)
  n <- nrow(Xtrain_norm)
  Kg <- kernel_matrix(Xquery_norm, Xtrain_norm[g_indices,, drop = FALSE], theta_g, global_kernel, isotropic)
  Kfull_g <- matrix(0, t, n)
  Kfull_g[, g_indices] <- Kg
  Kfull_l <- matrix(0, t, n)
  if (ncol(local_indices) > 0L) {
    for (i in seq_len(t)) {
      idx <- local_indices[i, ]
      idx <- idx[idx > 0]
      if (length(idx)) {
        Kfull_l[i, idx] <- kernel_matrix(Xquery_norm[i,, drop = FALSE], Xtrain_norm[idx,, drop = FALSE], theta_l, local_kernel, TRUE)
      }
    }
  }
  K <- Kfull_g + Kfull_l
  alpha0 <- get_prior(prior = prior, model = "DKP", r0 = r0, p0 = p0, Y = Y, K = K)
  alpha_n <- alpha0 + K %*% Y
  prob <- alpha_n / rowSums(alpha_n)
  list(alpha0 = alpha0, alpha_n = alpha_n, prob = prob, K = K, K_global = Kfull_g, K_local = Kfull_l)
}

test_that("twin_dkp_posterior_rcpp matches dense reference on small fixtures", {
  X <- matrix(c(0, 0.25, 0.5, 0.75, 1), ncol = 1)
  Y <- matrix(c(3, 1, 0,  1, 3, 1,  0, 2, 3,  2, 2, 1,  1, 1, 4), ncol = 3, byrow = TRUE)
  g_indices <- c(1L, 3L)
  local_indices <- matrix(c(2L, 4L, 2L, 4L, 4L, 5L, 2L, 5L, 2L, 4L), nrow = 5, byrow = TRUE)
  for (prior in c("noninformative", "fixed", "adaptive")) {
    p0 <- c(0, 0.4, 0.6)
    ref <- .twindkp_dense_reference(X, X, Y, g_indices, local_indices, 0.4, 0.3, "gaussian", "wendland", TRUE, prior, 2, p0)
    new <- twin_dkp_posterior_rcpp(X, X, Y, g_indices, local_indices, 0.4, 0.3, "gaussian", "wendland", TRUE, prior, 2, p0)
    expect_equal(new$alpha0, ref$alpha0, tolerance = 1e-8)
    expect_equal(new$alpha_n, ref$alpha_n, tolerance = 1e-8)
    expect_equal(new$prob, ref$prob, tolerance = 1e-8)
  }
})

test_that("fit_TwinDKP supports priors, anisotropy, l=0, diagnostics, and optimization smoke", {
  fx <- make_twindkp_data()
  fit_adapt <- fit_TwinDKP(fx$X, fx$Y, prior = "adaptive", theta_g = .4, theta_l = .3, g = 10, l = 5, twins = 1)
  expect_s3_class(fit_adapt, "TwinDKP")
  fit_non <- fit_TwinDKP(fx$X, fx$Y, prior = "noninformative", theta_g = .4, theta_l = .3, g = 10, l = 0, twins = 1)
  expect_s3_class(fit_non, "TwinDKP")
  expect_equal(fit_non$control$l, 0L)

  X2 <- cbind(seq(0, 1, length.out = 12), rep(c(0.2, 0.8), 6))
  Y2 <- fx$Y[seq_len(12), ]
  fit_aniso <- fit_TwinDKP(X2, Y2, prior = "fixed", r0 = 2, p0 = c(0, 0.5, 0.5), theta_g = c(.4, .5), theta_l = .3, g = 5, l = 2, twins = 1, isotropic = FALSE)
  expect_s3_class(fit_aniso, "TwinDKP")

  fit_diag <- fit_TwinDKP(fx$X, fx$Y, prior = "fixed", p0 = c(0, 0.5, 0.5), theta_g = .4, theta_l = .3, g = 10, l = 5, twins = 1, store_kernel = TRUE)
  expect_equal(dim(fit_diag$diagnostics$K), c(nrow(fx$X), nrow(fx$X)))
  expect_equal(dim(fit_diag$diagnostics$K_global), c(nrow(fx$X), nrow(fx$X)))
  expect_equal(dim(fit_diag$diagnostics$K_local), c(nrow(fx$X), nrow(fx$X)))
  fit_nodiag <- make_twindkp_model_1d()$model
  expect_null(fit_nodiag$diagnostics$K)
  expect_null(fit_nodiag$diagnostics$K_global)
  expect_null(fit_nodiag$diagnostics$K_local)

  smoke <- fit_TwinDKP(
    fx$X[1:12, , drop = FALSE],
    fx$Y[1:12, ],
    prior = "fixed",
    p0 = c(0, 0.5, 0.5),
    theta_g = NULL,
    theta_l = .3,
    g = 5,
    l = 2,
    twins = 1,
    n_multi_start = 1,
    n_threads = 1
  )
  expect_s3_class(smoke, "TwinDKP")
  expect_equal(smoke$control$n_multi_start, 1L)
  expect_equal(smoke$control$n_threads, 1L)
})

test_that("fit_TwinDKP validates fixed p0 consistently", {
  fx <- make_twindkp_data()
  expect_s3_class(fit_TwinDKP(fx$X, fx$Y, prior = "fixed", p0 = c(0, 0.5, 0.5), theta_g = .4, theta_l = .3, g = 10, l = 5), "TwinDKP")
  expect_error(fit_TwinDKP(fx$X, fx$Y, prior = "fixed", p0 = c(.5, .5), theta_g = .4, theta_l = .3, g = 10, l = 5), "p0")
  expect_error(fit_TwinDKP(fx$X, fx$Y, prior = "fixed", p0 = c(.5, -.1, .6), theta_g = .4, theta_l = .3, g = 10, l = 5), "p0")
  expect_error(fit_TwinDKP(fx$X, fx$Y, prior = "fixed", p0 = c(.5, .3, .3), theta_g = .4, theta_l = .3, g = 10, l = 5), "p0")
  expect_error(fit_TwinDKP(fx$X, fx$Y, prior = "fixed", p0 = c(.5, NA, .5), theta_g = .4, theta_l = .3, g = 10, l = 5), "p0")
  expect_error(fit_TwinDKP(fx$X, fx$Y, prior = "fixed", p0 = c(.5, NaN, .5), theta_g = .4, theta_l = .3, g = 10, l = 5), "p0")
  expect_error(fit_TwinDKP(fx$X, fx$Y, prior = "fixed", p0 = c(.5, Inf, .5), theta_g = .4, theta_l = .3, g = 10, l = 5), "p0")
})
