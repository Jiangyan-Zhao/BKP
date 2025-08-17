set.seed(123)
n <- 5
q <- 3
X <- matrix(runif(n * 2), ncol = 2)
m_vec <- sample(5:10, n, replace = TRUE)
Y <- t(sapply(1:n, function(i) rmultinom(1, size = m_vec[i], prob = rep(1/q, q))))
K <- matrix(runif(n * n), nrow = n)

test_that("noninformative prior returns all ones", {
  alpha0 <- get_prior_dkp(prior = "noninformative", Y = Y, K = K)
  expect_equal(dim(alpha0), c(n, q))
  expect_true(all(alpha0 == 1))
})

test_that("fixed prior returns correct values", {
  r0 <- 2
  p0 <- rep(1/q, q)
  alpha0 <- get_prior_dkp(prior = "fixed", r0 = r0, p0 = p0, K = K)
  expect_equal(dim(alpha0), c(n, q))
  expect_true(all(alpha0 == r0 * p0))
})

test_that("adaptive prior computes alpha0 correctly", {
  r0 <- 2
  alpha0 <- get_prior_dkp(prior = "adaptive", r0 = r0, Y = Y, K = K)
  expect_equal(dim(alpha0), c(n, q))
  expect_true(all(alpha0 >= 1e-3))
})

test_that("adaptive prior errors on missing inputs", {
  expect_error(get_prior_dkp(prior = "adaptive"), "Either Y or p0 must be provided to determine class dimension q.")
  expect_error(get_prior_dkp(prior = "adaptive", Y = Y), "Y and K must be provided for adaptive prior.")
})

test_that("fixed prior errors on invalid inputs", {
  r0 <- -1
  p0 <- rep(1/q, q)
  expect_error(get_prior_dkp(prior = "fixed", r0 = r0, p0 = p0, K = K), "r0 must be positive")

  expect_error(get_prior_dkp(prior = "fixed", r0 = 1, p0 = c(0.5, 0.5), Y = Y, K = K),
               "Length of p0 must match the number of classes")

  expect_error(get_prior_dkp(prior = "fixed", r0 = 1, p0 = c(-0.1, 0.5, 0.6), K = K),
               "p0 must be a valid probability vector")
})
