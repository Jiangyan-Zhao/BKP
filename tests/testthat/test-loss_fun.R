
# ------------------ BKP (Binary) Example ------------------
test_that("BKP losses return numeric values", {
  set.seed(123)
  n <- 10
  d <- 2
  Xnorm <- matrix(runif(n * d), ncol = d)
  m <- rep(10, n)
  y <- rbinom(n, size = m, prob = runif(n))
  gamma <- rep(0, d)

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
  gamma <- rep(0, d)

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
