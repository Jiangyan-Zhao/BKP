make_bkp_model_1d <- function() {
  set.seed(123)

  true_pi_fun <- function(x) {
    (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
  }

  n <- 30
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- matrix(runif(n, Xbounds[1], Xbounds[2]), ncol = 1)
  true_pi <- true_pi_fun(X)
  m <- sample(100, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = true_pi)

  model <- fit_BKP(
    X, y, m,
    Xbounds = Xbounds,
    prior = "noninformative"
  )

  list(
    X = X,
    y = y,
    m = m,
    Xbounds = Xbounds,
    model = model
  )
}


make_bkp_model_3d_classification <- function() {
  set.seed(100)

  X <- matrix(runif(36), ncol = 3)
  m <- rep(1, nrow(X))
  y <- rbinom(nrow(X), size = 1, prob = 0.5)

  model <- fit_BKP(
    X, y, m,
    prior = "fixed",
    r0 = 5,
    p0 = 0.5,
    theta = c(0.3, 0.4, 0.5),
    isotropic = FALSE
  )

  list(
    X = X,
    y = y,
    m = m,
    model = model
  )
}
