
set.seed(123)
n <- 5
X <- matrix(runif(n * 2), ncol = 2)
y <- rbinom(n, size = 5, prob = 0.6)
m <- rep(5, n)
K <- matrix(runif(n * n), nrow = n)

test_that("noninformative prior returns all ones", {
  prior_res <- get_prior(prior = "noninformative", K = K)
  expect_equal(prior_res$alpha0, rep(1, n))
  expect_equal(prior_res$beta0, rep(1, n))
})

test_that("fixed prior returns correct values", {
  r0 <- 2
  p0 <- 0.3
  prior_res <- get_prior(prior = "fixed", r0 = r0, p0 = p0, K = K)
  expect_equal(prior_res$alpha0, rep(r0 * p0, n))
  expect_equal(prior_res$beta0, rep(r0 * (1 - p0), n))
})

test_that("adaptive prior computes alpha0 and beta0 correctly", {
  r0 <- 2
  prior_res <- get_prior(prior = "adaptive", r0 = r0, y = y, m = m, K = K)

  expect_length(prior_res$alpha0, n)
  expect_length(prior_res$beta0, n)
  expect_true(all(prior_res$alpha0 >= 1e-3))
  expect_true(all(prior_res$beta0 >= 1e-3))
})

test_that("adaptive prior errors on missing inputs", {
  expect_error(get_prior(prior = "adaptive"), "y, m, and K must be provided")
  expect_error(get_prior(prior = "adaptive", y = y, m = m), "y, m, and K must be provided")
})

test_that("fixed prior errors on invalid r0 or p0", {
  expect_error(get_prior(prior = "fixed", r0 = -1, p0 = 0.5, K = K), "r0 must be positive")
  expect_error(get_prior(prior = "fixed", r0 = 1, p0 = 0, K = K), "p0 must be in \\(0, 1\\)")
  expect_error(get_prior(prior = "fixed", r0 = 1, p0 = 1, K = K), "p0 must be in \\(0, 1\\)")
})
