test_that("loss_fun_dkp computes valid loss values for brier and log_loss", {
  # Generate toy data
  set.seed(123)
  n = 10
  Xnorm = matrix(runif(2 * n), ncol = 2)
  m = rep(10, n)
  y = rbinom(n, size = m, prob = runif(n))
  Y = cbind(y, m-y)

  gamma <- rep(0, 2)

  # Brier score
  loss1 <- loss_fun_dkp(gamma = gamma, Xnorm = Xnorm, Y = Y,
                        prior = "adaptive",
                        loss = "brier",
                        kernel = "matern32")
  expect_type(loss1, "double")
  expect_gte(loss1, 0)

  # Log-loss
  loss2 <- loss_fun_dkp(gamma = gamma, Xnorm = Xnorm, Y = Y,
                        loss = "log_loss",
                        kernel = "gaussian")
  expect_type(loss2, "double")
})
