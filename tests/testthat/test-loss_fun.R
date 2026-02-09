
# ------------------ BKP (Binary) Example ------------------
test_that("BKP losses return numeric values", {
  set.seed(123)
  n <- 10
  d <- 2
  Xnorm <- matrix(runif(n * d), ncol = d)
  m <- rep(10, n)
  y <- rbinom(n, size = m, prob = runif(n))
  gamma <- 0

  loss_types <- c("brier", "log_loss")

  for (loss_type in loss_types) {
    l <- loss_fun(
      gamma = gamma, Xnorm = Xnorm,
      y = y, m = m,
      model = "BKP",
      prior = "noninformative",
      loss = loss_type,
      kernel = "gaussian"
    )
    expect_true(is.numeric(l))
    expect_equal(length(l), 1)
    expect_false(is.na(l))
    expect_false(is.infinite(l))
  }
})

# ------------------ DKP (Multi-class) Example ------------------
test_that("DKP losses return numeric values", {
  set.seed(123)
  n <- 10
  q <- 3
  d <- 2
  Xnorm <- matrix(runif(n * d), ncol = d)
  Y <- t(rmultinom(n, size = 10, prob = rep(1/q, q))) # n x q
  gamma <- 0

  loss_types <- c("brier", "log_loss")

  for (loss_type in loss_types) {
    l <- loss_fun(
      gamma = gamma, Xnorm = Xnorm,
      Y = Y,
      model = "DKP",
      prior = "noninformative",
      loss = loss_type,
      kernel = "gaussian"
    )
    expect_true(is.numeric(l))
    expect_equal(length(l), 1)
    expect_false(is.na(l))
    expect_false(is.infinite(l))
  }
})


test_that("loss_fun supports isotropic lengthscale", {
  set.seed(123)
  n <- 10
  d <- 2
  Xnorm <- matrix(runif(n * d), ncol = d)
  y <- rbinom(n, size = 10, prob = 0.5)
  m <- rep(10, n)

  l <- loss_fun(
    gamma = 0,
    Xnorm = Xnorm,
    y = y,
    m = m,
    model = "BKP",
    prior = "noninformative",
    loss = "brier",
    kernel = "gaussian",
    isotropic = TRUE
  )

  expect_true(is.numeric(l))
  expect_length(l, 1)
})

test_that("loss_fun additional validation and DKP anisotropic branch", {
  set.seed(1)
  Xnorm <- matrix(runif(8), ncol = 2)
  y <- c(1, 0, 1, 0)
  m <- rep(1, 4)
  Y <- cbind(y, 1 - y)

  expect_error(loss_fun(gamma = "a", Xnorm = Xnorm, y = y, m = m, model = "BKP"), "'gamma' must be a numeric vector.")
  expect_error(loss_fun(gamma = 0, Xnorm = c(1, 2), y = y, m = m, model = "BKP"), "'Xnorm' must be a numeric matrix with no NA.")
  Xnorm_bad <- Xnorm
  Xnorm_bad[1, 1] <- NA
  expect_error(loss_fun(gamma = 0, Xnorm = Xnorm_bad, y = y, m = m, model = "BKP"), "'Xnorm' must be a numeric matrix with no NA.")

  expect_error(loss_fun(gamma = 0, Xnorm = Xnorm, y = y, m = m, model = "BKP", r0 = 0), "'r0' must be a positive scalar.")
  expect_error(loss_fun(gamma = 0, Xnorm = Xnorm, y = y, m = m, model = "BKP", p0 = -1), "'p0' must be numeric and nonnegative.")
  expect_error(loss_fun(gamma = 0, Xnorm = Xnorm, y = y, m = m, model = "BKP", isotropic = c(TRUE, FALSE)), "'isotropic' must be a single logical value.")

  l_dkp <- loss_fun(gamma = c(-1, -1), Xnorm = Xnorm, Y = Y, model = "DKP", kernel = "matern52", isotropic = FALSE)
  expect_true(is.numeric(l_dkp) && length(l_dkp) == 1)
})
