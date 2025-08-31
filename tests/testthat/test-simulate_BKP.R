test_that("simulate.BKP returns expected structure and values", {
  # 0. Define true success probability function
  true_pi_fun <- function(x) {
    (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
  }

  # 1. Simulate binary outcome data
  set.seed(123)
  n <- 30
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- tgp::lhs(n = n, rect = Xbounds)
  true_pi <- true_pi_fun(X)
  m <- sample(100, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = true_pi)

  # 2. Fit BKP model
  model <- fit.BKP(X, y, m, Xbounds = Xbounds)

  # 3. Make predictions
  n_Xnew <- 10
  Xnew <- matrix(seq(-2, 2, length.out = n_Xnew), ncol = 1)

  # 4. Run simulate
  nsim <- 5
  sim_result <- simulate(model, Xnew = Xnew, nsim = nsim, threshold = 0.5)

  # 5. Check class and structure
  expect_type(sim_result, "list")
  expect_in(names(sim_result), c("samples", "mean", "class", "X", "Xnew", "threshold"))

  # sims: matrix [n_new x nsim]
  expect_true(is.matrix(sim_result$samples))
  expect_equal(dim(sim_result$samples), c(n_Xnew, nsim))

  # mean: numeric vector
  expect_true(is.numeric(sim_result$mean))
  expect_length(sim_result$mean, n_Xnew)

  # class: binary matrix
  expect_true(is.matrix(sim_result$class))
  expect_true(all(sim_result$class %in% c(0L, 1L)))

  # Each column of class corresponds to a simulation
  expect_equal(dim(sim_result$class), c(n_Xnew, nsim))
})
