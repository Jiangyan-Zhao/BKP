# Test file for the fit_TwinBKP function
#
# The tests intentionally use fixed theta_g in most cases to keep runtime small.
# Hyperparameter optimization is covered by one small smoke test.

test_that("fit_TwinBKP runs correctly with a documentation-style example", {
  set.seed(1)

  n <- 30
  X <- matrix(runif(n * 2), ncol = 2)
  m <- sample(10:30, n, replace = TRUE)
  p <- plogis(2 * (X[, 1] - 0.5))
  y <- rbinom(n, size = m, prob = p)

  expect_no_error({
    model <- fit_TwinBKP(X, y, m, g = 10, twins = 2, theta_g = 0.25)
  })

  expect_s3_class(model, "TwinBKP")
  expect_true(is.list(model))
  expect_true(all(model$global_indices >= 1L))
  expect_true(all(model$global_indices <= nrow(X)))
  expect_equal(length(model$global_indices), length(unique(model$global_indices)))
  expect_true(is.finite(model$theta_l))
  expect_true(model$theta_l > 0)
  expect_equal(model$theta_g, 0.25)
})


test_that("fit_TwinBKP handles input validation similarly to fit_BKP", {
  expect_error(fit_TwinBKP(), "Arguments 'X', 'y', and 'm' must be provided.")

  n <- 12
  d <- 2
  X_test <- matrix(runif(n * d), nrow = n)
  y_test <- rbinom(n, size = 20, prob = 0.5)
  m_test <- rep(20, n)

  expect_error(
    fit_TwinBKP(X = X_test, y = y_test),
    "Arguments 'X', 'y', and 'm' must be provided."
  )

  expect_error(
    fit_TwinBKP(X = as.character(X_test), y = y_test, m = m_test),
    "'X' must be a numeric matrix or data frame."
  )

  expect_error(
    fit_TwinBKP(X = X_test, y = as.character(y_test), m = m_test),
    "'y' must be numeric."
  )

  expect_error(
    fit_TwinBKP(X = X_test, y = y_test, m = as.character(m_test)),
    "'m' must be numeric."
  )

  expect_error(
    fit_TwinBKP(X = X_test, y = y_test[1:(n - 1)], m = m_test),
    "'y' must have the same number of rows as 'X'."
  )

  expect_error(
    fit_TwinBKP(X = X_test, y = y_test, m = m_test[1:(n - 1)]),
    "'m' must have the same number of rows as 'X'."
  )

  expect_error(
    fit_TwinBKP(X = X_test, y = rep(-1, n), m = m_test),
    "'y' must be nonnegative."
  )

  expect_error(
    fit_TwinBKP(X = X_test, y = y_test, m = rep(0, n)),
    "'m' must be strictly positive."
  )

  expect_error(
    fit_TwinBKP(X = X_test, y = m_test + 1, m = m_test),
    "Each element of 'y' must be less than or equal to corresponding element of 'm'."
  )

  X_na <- X_test
  X_na[1, 1] <- NA
  expect_error(
    fit_TwinBKP(X = X_na, y = y_test, m = m_test),
    "Missing values are not allowed in 'X', 'y', or 'm'."
  )

  X_inf <- X_test
  X_inf[1, 1] <- Inf
  expect_error(
    fit_TwinBKP(X = X_inf, y = y_test, m = m_test),
    "'X', 'y', and 'm' must contain only finite values."
  )

  expect_error(
    fit_TwinBKP(X_test, y_test, m_test, store_kernel = c(TRUE, FALSE)),
    "'store_kernel' must be a single logical value."
  )
})


test_that("fit_TwinBKP validates Xbounds and prior settings", {
  set.seed(2)

  n <- 12
  d <- 2
  X <- matrix(runif(n * d), nrow = n)
  y <- rbinom(n, size = 20, prob = 0.5)
  m <- rep(20, n)

  expect_error(
    fit_TwinBKP(X, y, m, Xbounds = 1),
    "'Xbounds' must be a numeric matrix."
  )

  expect_error(
    fit_TwinBKP(X, y, m, Xbounds = matrix(1, nrow = d, ncol = 3)),
    paste0("'Xbounds' must be a matrix with dimensions d x 2, where d = ", d, ".")
  )

  expect_error(
    fit_TwinBKP(X, y, m, Xbounds = matrix(c(0, 1, 0.5, 0.4), ncol = 2, byrow = TRUE)),
    "Each row of 'Xbounds' must satisfy lower < upper"
  )

  expect_error(
    fit_TwinBKP(X, y, m, r0 = 0),
    "'r0' must be a positive scalar"
  )

  expect_error(
    fit_TwinBKP(X, y, m, prior = "fixed", p0 = 1),
    "For fixed prior in BKP, 'p0' must be a scalar in \\(0, 1\\)."
  )

  X_out <- matrix(runif(n * d, min = -0.2, max = 1.2), nrow = n)
  expect_warning(
    expect_s3_class(
      fit_TwinBKP(X_out, y, m, theta_g = 0.25, g = 6, twins = 1),
      "TwinBKP"
    ),
    "Input X does not appear to be normalized to \\[0,1\\]"
  )
})


test_that("fit_TwinBKP validates hyperparameter and twinning controls", {
  set.seed(3)

  n <- 14
  d <- 3
  X <- matrix(runif(n * d), nrow = n)
  y <- rbinom(n, size = 15, prob = 0.5)
  m <- rep(15, n)

  expect_error(
    fit_TwinBKP(X, y, m, theta_g = "bad"),
    "'theta_g' must be numeric."
  )

  expect_error(
    fit_TwinBKP(X, y, m, theta_g = c(0.2, 0.3), isotropic = TRUE),
    "When isotropic=TRUE, 'theta_g' must be a scalar."
  )

  expect_error(
    fit_TwinBKP(X, y, m, theta_g = c(0.2, 0.3), isotropic = FALSE),
    "length 3"
  )

  expect_error(
    fit_TwinBKP(X, y, m, theta_g = -0.2),
    "'theta_g' must be strictly positive."
  )

  expect_error(
    fit_TwinBKP(X, y, m, theta_l = 0),
    "'theta_l' must be a positive scalar."
  )

  expect_error(
    fit_TwinBKP(X, y, m, isotropic = c(TRUE, FALSE)),
    "'isotropic' must be a single logical value."
  )

  expect_error(
    fit_TwinBKP(X, y, m, n_multi_start = 0),
    "'n_multi_start' must be a positive integer."
  )

  expect_error(
    fit_TwinBKP(X, y, m, n_threads = 0),
    "'n_threads' must be a positive integer."
  )

  expect_error(
    fit_TwinBKP(X, y, m, theta_g = 0.3, g = 6, twins = 0),
    "'twins' must be a positive integer."
  )

  expect_error(
    fit_TwinBKP(X, y, m, theta_g = 0.3, g = 6, twins = NA),
    "'twins' must be a positive integer."
  )

  expect_error(
    fit_TwinBKP(X, y, m, g = 1),
    "'g' must be an integer between 2 and n - 1."
  )

  expect_error(
    fit_TwinBKP(X, y, m, l = -1),
    "'l' must be a nonnegative integer."
  )

  expect_error(
    fit_TwinBKP(X, y, m, l = n),
    "'l' cannot exceed the number of non-global training points."
  )

})


test_that("fit_TwinBKP returns an object with expected structure and content", {
  set.seed(4)

  n <- 24
  d <- 2
  X <- matrix(runif(n * d), nrow = n)
  y <- rbinom(n, size = 12, prob = 0.45)
  m <- rep(12, n)

  model <- fit_TwinBKP(
    X = X, y = y, m = m,
    theta_g = 0.3, theta_l = 0.4,
    g = 8, twins = 2
  )

  expect_s3_class(model, "TwinBKP")
  expect_true(is.list(model))

  expected_names <- c(
    "theta_opt", "theta_g", "theta_l", "kernel", "global_kernel", "local_kernel",
    "isotropic", "loss", "loss_min", "ess", "ess_info",
    "X", "Xnorm", "Xbounds", "y", "m", "prior", "r0", "p0",
    "alpha0", "beta0", "alpha_n", "beta_n",
    "K", "K_global", "K_local", "twin_data", "twin_info",
    "global_indices", "local_indices", "g_target", "g", "r", "l", "twins", "u1",
    "leaf_size", "store_kernel", "complexity"
  )

  expect_true(all(expected_names %in% names(model)))

  expect_equal(model$theta_opt, model$theta_g)
  expect_equal(model$theta_g, 0.3)
  expect_equal(model$theta_l, 0.4)
  expect_equal(model$kernel, model$global_kernel)
  expect_equal(model$global_kernel, "gaussian")
  expect_equal(model$local_kernel, "wendland")
  expect_equal(model$loss, "brier")
  expect_equal(model$ess, "none")
  expect_equal(model$prior, "noninformative")
  expect_equal(model$twins, 2L)
  expect_true(is.integer(model$u1))
  expect_equal(length(model$u1), 2L)
  expect_true(is.numeric(model$r) || is.integer(model$r))
  expect_true(model$r >= 2L)
  expect_equal(model$leaf_size, 8L)
  expect_equal(model$X, X)
  expect_equal(dim(model$Xnorm), dim(X))
  expect_null(model$K)
  expect_null(model$K_global)
  expect_null(model$K_local)
  expect_false(model$store_kernel)
  expect_true(all(is.finite(model$alpha_n)))
  expect_true(all(is.finite(model$beta_n)))
  expect_equal(nrow(model$local_indices), nrow(model$X))
  expect_equal(ncol(model$local_indices), model$l)
  expect_false(any(model$local_indices %in% model$global_indices))
})


test_that("fit_TwinBKP uses fixed augmented Twinning data", {
  set.seed(123)

  n <- 20
  X <- matrix(runif(n * 2), ncol = 2)
  m <- sample(8:20, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = 0.4)

  model <- fit_TwinBKP(
    X, y, m,
    theta_g = 0.3,
    theta_l = 0.4,
    g = 6,
    twins = 2
  )

  expect_equal(ncol(model$twin_data), ncol(model$Xnorm) + 1L)
  expect_equal(
    model$twin_data[, ncol(model$twin_data)],
    as.numeric(y) / as.numeric(m)
  )
})


test_that("fit_TwinBKP stores dense diagnostic kernels only on request", {
  set.seed(44)

  n <- 24
  X <- matrix(runif(n * 2), ncol = 2)
  m <- sample(8:15, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = 0.45)

  model <- fit_TwinBKP(
    X = X, y = y, m = m,
    theta_g = 0.3, theta_l = 0.4,
    g = 8, twins = 2, store_kernel = TRUE,
    ess = "none"
  )

  expect_true(is.matrix(model$K))
  expect_true(is.matrix(model$K_global))
  expect_true(is.matrix(model$K_local))
  expect_equal(model$K, model$K_global + model$K_local)
  expect_equal(model$alpha_n, model$alpha0 + as.vector(model$K %*% model$y))
  expect_equal(model$beta_n, model$beta0 + as.vector(model$K %*% (model$m - model$y)))
})


test_that("fit_TwinBKP computes theta_l as the empirical covering radius", {
  set.seed(5)

  n <- 30
  X <- matrix(runif(n * 2), ncol = 2)
  y <- rbinom(n, size = 10, prob = 0.5)
  m <- rep(10, n)

  model <- fit_TwinBKP(
    X = X, y = y, m = m,
    theta_g = 0.25, g = 10, twins = 2
  )

  G <- model$global_indices
  theta_manual <- max(vapply(seq_len(nrow(model$Xnorm)), function(i) {
    dif <- sweep(model$Xnorm[G, , drop = FALSE], 2, model$Xnorm[i, ], "-")
    min(sqrt(rowSums(dif^2)))
  }, numeric(1)))

  expect_equal(model$theta_l, theta_manual, tolerance = 1e-10)
})


test_that("fit_TwinBKP posterior matches the combined global-local kernel update", {
  set.seed(6)

  n <- 22
  X <- matrix(runif(n * 2), ncol = 2)
  m <- sample(8:15, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = 0.4)

  model <- fit_TwinBKP(
    X = X, y = y, m = m,
    theta_g = 0.25, theta_l = 0.35,
    g = 8, twins = 2,
    ess = "none",
    store_kernel = TRUE
  )

  prior_par <- get_prior(
    prior = model$prior,
    model = "BKP",
    r0 = model$r0,
    p0 = model$p0,
    y = model$y,
    m = model$m,
    K = model$K
  )

  expect_equal(model$alpha0, prior_par$alpha0)
  expect_equal(model$beta0, prior_par$beta0)
  expect_equal(model$alpha_n, prior_par$alpha0 + as.vector(model$K %*% model$y))
  expect_equal(model$beta_n, prior_par$beta0 + as.vector(model$K %*% (model$m - model$y)))
  expect_equal(model$ess_info$scale, rep(1, n))
})


test_that("fit_TwinBKP supports fixed prior, adaptive prior, and anisotropic theta_g", {
  set.seed(7)

  n <- 24
  X <- matrix(runif(n * 2), ncol = 2)
  m <- sample(8:20, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = 0.45)

  fixed_model <- fit_TwinBKP(
    X, y, m,
    prior = "fixed", r0 = 0.5, p0 = 0.5,
    theta_g = 0.3, theta_l = 0.4,
    g = 8, twins = 2
  )

  adaptive_model <- fit_TwinBKP(
    X, y, m,
    prior = "adaptive", r0 = mean(m),
    theta_g = c(0.3, 0.4), theta_l = 0.4,
    isotropic = FALSE,
    g = 8, twins = 2
  )

  expect_s3_class(fixed_model, "TwinBKP")
  expect_s3_class(adaptive_model, "TwinBKP")
  expect_equal(length(adaptive_model$theta_g), 2)
  expect_equal(adaptive_model$isotropic, FALSE)
})


test_that("fit_TwinBKP handles Shepard ESS without NA values", {
  set.seed(8)

  X <- matrix(c(
    0.05, 0.15,
    0.20, 0.70,
    0.35, 0.30,
    0.50, 0.85,
    0.65, 0.20,
    0.80, 0.60,
    0.95, 0.40,
    0.12, 0.92
  ), ncol = 2, byrow = TRUE)
  m <- c(10, 14, 18, 12, 20, 16, 22, 15)
  y <- c(3, 6, 8, 7, 11, 9, 15, 5)

  model <- fit_TwinBKP(
    X, y, m,
    theta_g = 0.35,
    theta_l = 0.45,
    g = 4,
    twins = 2,
    ess = "shepard"
  )

  expect_s3_class(model, "TwinBKP")
  expect_false(anyNA(model$alpha_n))
  expect_false(anyNA(model$beta_n))
  expect_false(anyNA(model$ess_info$scale))
  expect_true(all(is.finite(model$alpha_n)))
  expect_true(all(is.finite(model$beta_n)))
  expect_true(all(is.finite(model$ess_info$scale)))
})


test_that("fit_TwinBKP rejects duplicated input locations for Shepard ESS", {
  X <- matrix(c(0.10, 0.20, 0.10, 0.20, 0.75, 0.90, 0.40, 0.55),
              ncol = 2, byrow = TRUE)
  m <- c(10, 20, 15, 12)
  y <- c(4, 9, 7, 6)

  expect_error(
    fit_TwinBKP(X, y, m, theta_g = 0.3, ess = "shepard", g = 2, twins = 1),
    "requires unique input locations"
  )
})


test_that("fit_TwinBKP is reproducible with a fixed random seed", {
  set.seed(9)

  n <- 26
  X <- matrix(runif(n * 2), ncol = 2)
  m <- sample(10:20, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = 0.4)

  set.seed(123)
  model1 <- fit_TwinBKP(
    X, y, m,
    theta_g = 0.25,
    g = 8,
    twins = 3
  )

  set.seed(123)
  model2 <- fit_TwinBKP(
    X, y, m,
    theta_g = 0.25,
    g = 8,
    twins = 3
  )

  expect_equal(model1$u1, model2$u1)
  expect_equal(model1$global_indices, model2$global_indices)
  expect_equal(model1$theta_l, model2$theta_l)
  expect_equal(model1$K, model2$K)
})


test_that("fit_TwinBKP can optimize theta_g on a small global subset", {
  set.seed(10)

  n <- 18
  X <- matrix(runif(n * 2), ncol = 2)
  m <- sample(8:15, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = 0.5)

  model <- fit_TwinBKP(
    X, y, m,
    g = 6,
    twins = 1,
    n_multi_start = 1,
    n_threads = 1
  )

  expect_s3_class(model, "TwinBKP")
  expect_true(all(is.finite(model$theta_g)))
  expect_true(all(model$theta_g > 0))
  expect_true(is.finite(model$loss_min))
})
