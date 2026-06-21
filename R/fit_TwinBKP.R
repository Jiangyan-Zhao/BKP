fit_TwinBKP <- function(
    X, y, m, Xbounds = NULL,
    prior = c("noninformative", "fixed", "adaptive"),
    r0 = 2, p0 = mean(y/m),
    global_kernel = c("gaussian", "matern52", "matern32", "wendland"),
    local_kernel = c("wendland"),
    loss = c("brier", "log_loss"),
    n_multi_start = NULL,
    theta_g = NULL,
    theta_l = NULL,
    isotropic = TRUE,
    n_threads = 1,
    ess = c("none", "shepard"),
    g = NULL,
    r = NULL,
    l = NULL,
    runs = 10,
    u1 = NULL,
    leaf_size = 8,
    response_weight = 1,
    include_m_in_twin = FALSE,
    size_weight = 0
) {

  prior <- match.arg(prior)
  global_kernel <- match.arg(global_kernel)
  local_kernel <- match.arg(local_kernel)
  loss <- match.arg(loss)
  ess <- match.arg(ess)

  X <- as.matrix(X)
  y <- as.numeric(y)
  m <- as.numeric(m)

  n <- nrow(X); d <- ncol(X)

  ## ---------------- normalize ----------------
  if (is.null(Xbounds)) {
    Xbounds <- cbind(rep(0, d), rep(1, d))
  }

  Xnorm <- sweep(X, 2, Xbounds[,1], "-")
  Xnorm <- sweep(Xnorm, 2, Xbounds[,2] - Xbounds[,1], "/")

  ## ---------------- twin_data ----------------
  p_hat <- y / m
  twin_data <- cbind(Xnorm, response_weight * p_hat)

  ## ---------------- g / r ----------------
  if (is.null(g)) {
    g <- min(floor(n / 2), 50 * d, max(floor(sqrt(n)), 10 * d))
    g <- max(2L, min(n - 1L, g))
  }

  if (is.null(r)) {
    r <- ceiling(n / g)
  }

  r <- max(2L, as.integer(r))

  ## ---------------- starts ----------------
  if (is.null(u1)) {
    center <- colMeans(twin_data)
    d2 <- rowSums((twin_data - center)^2)
    u1 <- c(which.max(d2),
            sample.int(n, runs - 1, replace = TRUE))
  }

  ## ---------------- C++ twinning ----------------
  twin_info <- twin_select_global_rcpp(
    twin_data, Xnorm,
    r = r, runs = runs,
    u1 = as.integer(u1),
    leaf_size = leaf_size
  )

  G <- as.integer(twin_info$g_indices)

  ## ---------------- local neighbors ----------------
  nonG <- setdiff(seq_len(n), G)
  l <- if (is.null(l)) min(length(nonG), max(25L, 3 * d)) else l

  Lmat <- matrix(NA_integer_, n, l)
  for (i in seq_len(n)) {
    if (l > 0 && length(nonG) > 0) {
      d2 <- rowSums((Xnorm[nonG,,drop=FALSE] - Xnorm[i,])^2)
      Lmat[i,] <- nonG[order(d2)[seq_len(l)]]
    }
  }

  ## ---------------- theta_l ----------------
  theta_l <- as.numeric(twin_info$theta_l)
  if (!is.finite(theta_l) || theta_l <= 0)
    stop("Invalid theta_l")

  ## ---------------- global theta ----------------
  Xg <- Xnorm[G,,drop=FALSE]
  yg <- y[G]; mg <- m[G]

  if (is.null(theta_g)) {
    fit_g <- fit_BKP(
      X = X[g,,drop=FALSE],
      y = yg,
      m = mg,
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
      ess = ess
    )
    theta_g <- fit_g$theta_opt
  }

  ## ---------------- kernel construction ----------------
  Kg <- matrix(0, n, n)
  Kg[,G] <- kernel_matrix(Xnorm, Xg,
                          theta = theta_g,
                          kernel = global_kernel,
                          isotropic = isotropic)

  Kl <- matrix(0, n, n)
  if (l > 0) {
    for (i in seq_len(n)) {
      idx <- Lmat[i,]
      idx <- idx[!is.na(idx)]
      if (length(idx) > 0) {
        Kl[i, idx] <- kernel_matrix(
          Xnorm[i,,drop=FALSE],
          Xnorm[idx,,drop=FALSE],
          theta = theta_l,
          kernel = "wendland",
          isotropic = TRUE
        )
      }
    }
  }

  K <- Kg + Kl

  ## ---------------- posterior ----------------
  prior_par <- get_prior(prior, "BKP",
                         r0 = r0, p0 = p0,
                         y = y, m = m, K = K)

  alpha0 <- prior_par$alpha0
  beta0 <- prior_par$beta0

  data_success <- as.vector(K %*% y)
  data_failure <- as.vector(K %*% (m - y))

  if (identical(ess, "shepard")) {
    ess_info <- .bkp_ess_calibration(Xnorm, Xnorm, m, K)
    data_success <- ess_info$scale * data_success
    data_failure <- ess_info$scale * data_failure
  } else {
    ess_info <- .bkp_ess_none_info(K, m)
  }

  structure(list(
    theta_g = theta_g,
    theta_l = theta_l,
    global_indices = G,
    local_indices = Lmat,
    K = K,
    alpha0 = alpha0,
    beta0 = beta0,
    alpha_n = alpha0 + data_success,
    beta_n = beta0 + data_failure,
    ess_info = ess_info,
    X = X,
    Xnorm = Xnorm
  ), class = "TwinBKP")
}
