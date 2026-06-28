#' Fit a Twin Dirichlet Kernel Process Model
#'
#' Fits a scalable global-local approximation to [fit_DKP()] for categorical or
#' multinomial count responses.
#'
#' @inheritParams fit_DKP
#' @param theta_g Optional positive global kernel lengthscale(s).
#' @param theta_l Optional positive local kernel range.
#' @param global_kernel,local_kernel Kernel functions for global and local components.
#' @param g Target global subset size.
#' @param l Number of local non-global neighbours.
#' @param n_multi_start Optional number of multi-start initial points for global DKP hyperparameter tuning.
#' @param n_threads Number of threads used for global DKP hyperparameter tuning.
#' @param twins Number of Twinning runs.
#' @param store_kernel Store dense diagnostic kernels.
#'
#' @return A `TwinDKP` object.
#'
#' @examples
#' set.seed(2026)
#' X <- matrix(seq(0, 1, length.out = 20), ncol = 1)
#' P <- cbind(1 - X[, 1], rep(0.3, 20), X[, 1] + 0.1)
#' P <- P / rowSums(P)
#' Y <- t(vapply(seq_len(20), function(i) as.numeric(rmultinom(1, 8, P[i, ])), numeric(3)))
#' fit_TwinDKP(X, Y, prior = "fixed", p0 = rep(1/3, 3), theta_g = 0.4,
#'             theta_l = 0.3, g = 6, l = 4, twins = 1)
#'
#' @export
fit_TwinDKP <- function(
    X, Y, Xbounds = NULL,
    prior = c("noninformative", "fixed", "adaptive"),
    r0 = 2, p0 = NULL,
    global_kernel = c("gaussian", "matern52", "matern32", "wendland"),
    local_kernel = c("wendland", "gaussian", "matern52", "matern32"),
    loss = c("brier", "log_loss"),
    n_multi_start = NULL,
    isotropic = TRUE,
    n_threads = 1,
    theta_g = NULL,
    theta_l = NULL,
    g = NULL,
    l = NULL,
    twins = 5,
    store_kernel = FALSE
) {
  # ---- Argument checking ----
  if (missing(X) || missing(Y)) {
    stop("Arguments 'X' and 'Y' must be provided.")
  }
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("'X' must be a numeric matrix or data frame.")
  }
  if (!is.numeric(as.matrix(X))) {
    stop("'X' must contain numeric values only.")
  }
  if (!is.matrix(Y) && !is.data.frame(Y)) {
    stop("'Y' must be a numeric matrix or data frame.")
  }
  if (!is.numeric(as.matrix(Y))) {
    stop("'Y' must contain numeric values only.")
  }

  X <- as.matrix(X)
  Y <- as.matrix(Y)

  d <- ncol(X)
  n <- nrow(X)
  q <- ncol(Y)

  if (n < 3L) {
    stop("TwinDKP requires at least three observations.")
  }
  if (nrow(Y) != n) {
    stop("Number of rows in 'Y' must match number of rows in 'X'.")
  }
  if (q < 2L) {
    stop("'Y' must have at least two columns (multinomial outcomes).")
  }
  if (any(Y < 0)) {
    stop("'Y' must be nonnegative counts or frequencies.")
  }
  if (any(rowSums(Y) <= 0)) {
    stop("Each row of 'Y' must have a positive row sum.")
  }
  if (anyNA(X) || anyNA(Y)) {
    stop("Missing values are not allowed in 'X' or 'Y'.")
  }
  if (any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("'X' and 'Y' must contain only finite values.")
  }
  # ---- Prior, kernels, and loss ----
  prior <- match.arg(prior)
  global_kernel <- match.arg(global_kernel)
  local_kernel <- match.arg(local_kernel)
  loss <- match.arg(loss)

  if (is.null(p0)) {
    p0 <- colMeans(sweep(Y, 1, rowSums(Y), "/"))
  }

  if (!is.numeric(r0) || length(r0) != 1 || r0 <= 0 || is.na(r0) || !is.finite(r0)) {
    stop("'r0' must be a positive scalar.")
  }
  if (prior == "fixed" &&
      (is.null(p0) || !is.numeric(p0) || length(p0) != q || anyNA(p0) ||
       any(!is.finite(p0)) || any(p0 < 0) || abs(sum(p0) - 1) > 1e-10)) {
    stop("For fixed prior in DKP, 'p0' must be a nonnegative finite numeric vector of length equal to the number of classes and sum to 1.")
  }
  if (!is.logical(isotropic) || length(isotropic) != 1) {
    stop("'isotropic' must be a single logical value.")
  }

  if (!is.null(n_multi_start)) {
    if (!is.numeric(n_multi_start) || length(n_multi_start) != 1 ||
        is.na(n_multi_start) || !is.finite(n_multi_start) ||
        n_multi_start <= 0 || n_multi_start != floor(n_multi_start)) {
      stop("'n_multi_start' must be a positive integer.")
    }
    n_multi_start <- as.integer(n_multi_start)
  }

  if (!is.numeric(n_threads) || length(n_threads) != 1 ||
      is.na(n_threads) || !is.finite(n_threads) || n_threads <= 0 ||
      n_threads != floor(n_threads)) {
    stop("'n_threads' must be a positive integer.")
  }
  if (!is.null(theta_g)) {
    if (!is.numeric(theta_g)) {
      stop("'theta_g' must be numeric.")
    }
    if (isTRUE(isotropic) && length(theta_g) != 1) {
      stop("When isotropic=TRUE, 'theta_g' must be a scalar.")
    }
    if (!isTRUE(isotropic) && !(length(theta_g) == 1 || length(theta_g) == d)) {
      stop(paste0("When isotropic=FALSE, 'theta_g' must be either a scalar or a vector of length ", d, "."))
    }
    if (!isTRUE(isotropic) && length(theta_g) == 1) {
      theta_g <- rep(theta_g, d)
    }
    if (anyNA(theta_g) || any(!is.finite(theta_g)) || any(theta_g <= 0)) {
      stop("'theta_g' must be strictly positive.")
    }
  }
  if (!is.null(theta_l) &&
      (!is.numeric(theta_l) || length(theta_l) != 1 || is.na(theta_l) ||
       !is.finite(theta_l) || theta_l <= 0)) {
    stop("'theta_l' must be a positive scalar.")
  }
  if (!is.numeric(twins) || length(twins) != 1 || is.na(twins) ||
      !is.finite(twins) || twins <= 0 || twins != floor(twins)) {
    stop("'twins' must be a positive integer.")
  }
  if (!is.logical(store_kernel) || length(store_kernel) != 1) {
    stop("'store_kernel' must be a single logical value.")
  }

  # ---- Normalize input X to [0,1]^d ----
  if (is.null(Xbounds)) {
    xmin <- min(X)
    xmax <- max(X)
    if (xmin < 0 || xmax > 1) {
      warning(sprintf(
        "Input X does not appear to be normalized to [0,1]. Current range: [%.3f, %.3f].\nPlease normalize X or specify Xbounds explicitly; otherwise the model may produce incorrect results.",
        xmin, xmax
      ))
    }
    Xbounds <- cbind(rep(0, d), rep(1, d))
  } else {
    if (!is.matrix(Xbounds)) {
      stop("'Xbounds' must be a numeric matrix.")
    }
    if (!is.numeric(Xbounds)) {
      stop("'Xbounds' must contain numeric values.")
    }
    if (!all(dim(Xbounds) == c(d, 2))) {
      stop(paste0("'Xbounds' must be a matrix with dimensions d x 2, where d = ", d, "."))
    }
    if (anyNA(Xbounds) || any(!is.finite(Xbounds))) {
      stop("'Xbounds' must contain only finite values.")
    }
    if (any(Xbounds[, 2] <= Xbounds[, 1])) {
      stop("Each row of 'Xbounds' must satisfy lower < upper.")
    }
  }

  Xnorm <- sweep(X, 2, Xbounds[, 1], "-")
  Xnorm <- sweep(Xnorm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

  # ---- Twinning global subset selection ----
  if (is.null(g)) {
    g <- as.integer(max(2L, min(n - 1L, 50L * d, max(floor(sqrt(n)), 10L * d))))
  } else {
    if (!is.numeric(g) || length(g) != 1 || g != floor(g) ||
        is.na(g) || !is.finite(g) || g < 2 || g >= n) {
      stop("'g' must be an integer between 2 and n - 1.")
    }
    g <- as.integer(g)
  }

  r <- as.integer(max(2L, ceiling(n / g)))
  twins <- as.integer(twins)
  leaf_size <- 8L
  P_hat <- sweep(Y, 1, rowSums(Y), "/")
  twin_data <- cbind(Xnorm, P_hat)
  storage.mode(twin_data) <- "double"
  center <- colMeans(twin_data)
  first <- which.max(rowSums(sweep(twin_data, 2, center, "-")^2))
  u1 <- if (twins <= 1L) {
    as.integer(first)
  } else {
    as.integer(c(
      first,
      sample(
        setdiff(seq_len(n), first),
        twins - 1L,
        replace = (twins - 1L) > (n - 1L)
      )
    ))
  }

  twin_info <- twin_select_global_rcpp(
    twin_data = twin_data,
    Xnorm = Xnorm,
    r = r,
    runs = twins,
    u1 = u1,
    leaf_size = leaf_size
  )
  g_indices <- as.integer(twin_info$g_indices)
  g_actual <- length(g_indices)

  if (g_actual < 2L) {
    stop("The Twinning step selected fewer than two global points; try increasing 'g'.")
  }
  if (is.null(theta_l)) {
    theta_l <- as.numeric(twin_info$theta_l)
  }

  non_global_n <- n - g_actual
  if (is.null(l)) {
    l <- min(non_global_n, max(25L, 3L * d))
  } else {
    if (!is.numeric(l) || length(l) != 1 || l != floor(l) ||
        is.na(l) || !is.finite(l) || l < 0) {
      stop("'l' must be a nonnegative integer.")
    }
    l <- as.integer(l)
    if (l > non_global_n) {
      stop("'l' cannot exceed the number of non-global training points.")
    }
  }

  # ---- Global lengthscale tuning ----
  if (is.null(theta_g)) {
    global_fit <- fit_DKP(
      X[g_indices, , drop = FALSE],
      Y[g_indices, , drop = FALSE],
      Xbounds = Xbounds,
      prior = prior,
      r0 = r0,
      p0 = p0,
      kernel = global_kernel,
      loss = loss,
      n_multi_start = n_multi_start,
      theta = NULL,
      isotropic = isotropic,
      n_threads = n_threads,
      ess = "none"
    )
    theta_g <- as.numeric(global_fit$theta_opt)
    loss_min <- as.numeric(global_fit$loss_min)
  } else {
    theta_g <- as.numeric(theta_g)
    loss_min <- loss_fun(
      gamma = log10(theta_g),
      Xnorm = Xnorm[g_indices, , drop = FALSE],
      Y = Y[g_indices, , drop = FALSE],
      prior = prior,
      r0 = r0,
      p0 = p0,
      model = "DKP",
      loss = loss,
      kernel = global_kernel,
      isotropic = isotropic,
      ess = "none"
    )
  }

  local_indices <- twin_local_indices_rcpp(
    Xtrain_norm = Xnorm,
    Xquery_norm = Xnorm,
    g_indices = g_indices,
    l = l,
    leaf_size = leaf_size
  )

  # ---- Compute posterior parameters ----
  posterior <- twin_dkp_posterior_rcpp(
    Xquery_norm = Xnorm,
    Xtrain_norm = Xnorm,
    Y = Y,
    g_indices = as.integer(g_indices),
    local_indices = local_indices,
    theta_g = as.numeric(theta_g),
    theta_l = as.numeric(theta_l),
    global_kernel = global_kernel,
    local_kernel = local_kernel,
    isotropic = isTRUE(isotropic),
    prior = prior,
    r0 = r0,
    p0 = as.numeric(p0),
    store_kernel = store_kernel
  )

  # ---- Construct and return fitted model ----
  out <- list(
    theta_opt = theta_g,
    theta_g = theta_g,
    theta_l = theta_l,
    kernel = global_kernel,
    global_kernel = global_kernel,
    local_kernel = local_kernel,
    isotropic = isotropic,
    loss = loss,
    loss_min = loss_min,
    X = X,
    Xnorm = Xnorm,
    Xbounds = Xbounds,
    Y = Y,
    prior = prior,
    r0 = r0,
    p0 = p0,
    alpha0 = posterior$alpha0,
    alpha_n = posterior$alpha_n,
    prob = posterior$prob,
    global_indices = g_indices,
    control = list(
      g_target = g,
      g = g_actual,
      l = l,
      r = r,
      twins = twins,
      u1 = u1,
      leaf_size = leaf_size,
      n_multi_start = n_multi_start,
      n_threads = n_threads,
      store_kernel = store_kernel
    ),
    diagnostics = list(
      twin_info = twin_info,
      K = posterior$K,
      K_global = posterior$K_global,
      K_local = posterior$K_local
    ),
    ess = "none",
    ess_info = NULL
  )
  class(out) <- "TwinDKP"
  out
}
