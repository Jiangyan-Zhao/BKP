
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
