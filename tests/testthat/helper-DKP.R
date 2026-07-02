make_dkp_model_1d <- function() {
  set.seed(123)

  true_pi_fun <- function(X) {
    p1 <- 1 / (1 + exp(-3 * X))
    p2 <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
    matrix(c(p1 / 2, p2 / 2, 1 - (p1 + p2) / 2), nrow = length(p1))
  }

  n <- 30
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- matrix(runif(n, Xbounds[1], Xbounds[2]), ncol = 1)
  true_pi <- true_pi_fun(X)
  m <- sample(150, n, replace = TRUE)

  Y <- t(sapply(seq_len(n), function(i) {
    rmultinom(1, size = m[i], prob = true_pi[i, ])
  }))

  model <- fit_DKP(
    X, Y,
    Xbounds = Xbounds,
    prior = "noninformative",
    theta = 0.3
  )

  list(
    X = X,
    Y = Y,
    m = m,
    Xbounds = Xbounds,
    model = model
  )
}


make_dkp_model_3d_classification <- function() {
  set.seed(100)

  X <- matrix(runif(45), ncol = 3)
  Y <- matrix(0, nrow(X), 4)

  cls <- sample(1:4, nrow(X), replace = TRUE)
  Y[cbind(seq_len(nrow(X)), cls)] <- 1

  model <- fit_DKP(
    X, Y,
    prior = "fixed",
    r0 = 6,
    p0 = rep(0.25, 4),
    theta = c(0.35, 0.4, 0.45),
    isotropic = FALSE
  )

  list(
    X = X,
    Y = Y,
    model = model
  )
}
