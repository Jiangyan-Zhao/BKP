test_that("loss_fun computes valid loss values for brier and log_loss", {
  set.seed(123)
  n = 10
  Xnorm = matrix(runif(2 * n), ncol = 2)
  m = rep(10, n)
  y = rbinom(n, size = m, prob = runif(n))

  gamma <- rep(0, 2)

  # Brier score
  loss1 <- loss_fun(gamma = gamma, Xnorm = Xnorm, y = y, m = m)

  expect_type(loss1, "double")
  expect_gte(loss1, 0)

  # Log-loss
  loss2 <- loss_fun(gamma = gamma, Xnorm = Xnorm, y = y, m = m,
                    prior = "adaptive",
                    loss = "log_loss",
                    kernel = "matern52")
  expect_type(loss2, "double")
})
