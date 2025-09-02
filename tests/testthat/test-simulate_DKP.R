test_that("simulate.DKP returns expected structure and dimensions", {
  # 0. Define true class probability function (3-class)
  true_pi_fun <- function(X) {
    p <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
    return(matrix(c(p/2, p/2, 1 - p), nrow = length(p)))
  }

  # 1. Simulate training data
  n <- 30
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- tgp::lhs(n = n, rect = Xbounds)
  true_pi <- true_pi_fun(X)
  m <- rep(1, n)
  Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
  q <- ncol(Y)  # neumber of class

  # 2. Fit DKP model
  model <- fit_DKP(X, Y, Xbounds = Xbounds)

  # 3. New prediction locations
  n_Xnew <- 10
  Xnew <- matrix(seq(-2, 2, length.out = n_Xnew), ncol = 1)

  # 4. Simulate from the DKP posterior
  nsim <- 5
  sim_result <- simulate(model, Xnew = Xnew, nsim = nsim)

  # 5. Check output structure
  expect_type(sim_result, "list")
  expect_named(sim_result, c('samples', 'mean', 'class', 'X', 'Xnew'))

  # sims: array of dimension [nsim × q × n_new]
  expect_true(is.array(sim_result$samples))
  expect_equal(dim(sim_result$samples), c(n_Xnew, q, nsim))

  # mean: matrix of dimension [n_new × q]
  expect_true(is.matrix(sim_result$mean))
  expect_in(names(sim_result), c('samples', 'mean', 'class', 'X', 'Xnew'))

  # class: binary matrix
  expect_true(is.matrix(sim_result$class))
  expect_true(all(sim_result$class %in% 1:q))

  # Each column of class corresponds to a simulation
  expect_equal(dim(sim_result$class), c(n_Xnew, nsim))
})
