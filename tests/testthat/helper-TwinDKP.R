make_twindkp_data <- function() {
  set.seed(2026)
  n <- 40
  X <- matrix(seq(0, 1, length.out = n), ncol = 1)
  eta1 <- 1.5 - 2 * X[, 1]
  eta2 <- 0.2 + X[, 1]
  eta3 <- -0.5 + 1.5 * X[, 1]
  E <- cbind(eta1, eta2, eta3)
  P <- exp(E)
  P <- P / rowSums(P)
  m <- rep(12, n)
  Y <- t(vapply(seq_len(n), function(i) as.numeric(rmultinom(1, size = m[i], prob = P[i, ])), numeric(3)))
  list(X = X, Y = Y)
}
make_twindkp_model_1d <- function() {
  dat <- make_twindkp_data()
  model <- fit_TwinDKP(dat$X, dat$Y, prior = "fixed", r0 = 2, p0 = rep(1/3, 3), theta_g = 0.4, theta_l = 0.3, g = 10, l = 5, twins = 1, global_kernel = "gaussian", local_kernel = "wendland")
  c(dat, list(model = model))
}
