test_that("simulate.BKP returns expected structure and values", {
  skip_on_cran()

  # Generate toy data
  set.seed(123)
  n <- 20
  Xbounds <- matrix(c(-2, 2), nrow = 1)
  X <- tgp::lhs(n = n, rect = Xbounds)
  true_pi <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
  m <- sample(80:120, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = true_pi)

  # Fit model
  model <- fit.BKP(X, y, m, Xbounds = Xbounds)

  # Define new input
  Xnew <- matrix(seq(-2, 2, length.out = 10), ncol = 1)

  # Run simulate
  sim_result <- simulate(model, Xnew, n_sim = 5, threshold = 0.5)

  # Check class and structure
  expect_type(sim_result, "list")
  expect_named(sim_result, c("sims", "mean", "class"))

  # sims: matrix [n_new x n_sim]
  expect_true(is.matrix(sim_result$sims))
  expect_equal(dim(sim_result$sims), c(nrow(Xnew), 5))

  # mean: numeric vector
  expect_true(is.numeric(sim_result$mean))
  expect_length(sim_result$mean, nrow(Xnew))

  # class: binary matrix
  expect_true(is.matrix(sim_result$class))
  expect_equal(range(sim_result$class), c(0, 1))

  # Each column of class corresponds to a simulation
  expect_equal(dim(sim_result$class), c(nrow(Xnew), 5))
})
