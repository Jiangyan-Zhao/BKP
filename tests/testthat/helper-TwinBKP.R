make_twinbkp_model_1d <- function() {
  set.seed(123)
  n <- 30
  X <- matrix(seq(0, 1, length.out = n), ncol = 1)
  p <- plogis(8 * (X[, 1] - 0.5))
  m <- rep(20, n)
  y <- rbinom(n, size = m, prob = p)

  model <- fit_TwinBKP(
    X, y, m,
    Xbounds = matrix(c(0, 1), nrow = 1),
    theta_g = 0.25,
    theta_l = 0.30,
    g = 10,
    l = 5,
    twins = 1
  )

  list(X = X, y = y, m = m, model = model)
}

make_twinbkp_model_2d <- function() {
  set.seed(123)
  n <- 35
  X <- matrix(runif(2 * n), ncol = 2)
  p <- plogis(5 * (X[, 1] - X[, 2]))
  m <- rep(15, n)
  y <- rbinom(n, size = m, prob = p)

  model <- fit_TwinBKP(
    X, y, m,
    Xbounds = matrix(c(0, 0, 1, 1), nrow = 2),
    theta_g = 0.35,
    theta_l = 0.30,
    g = 10,
    l = 5,
    twins = 1
  )

  list(X = X, y = y, m = m, model = model)
}
