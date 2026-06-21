#' @name fit_TwinBKP
#'
#' @title Fit a Twin Beta Kernel Process Model
#'
#' @description Fits a Twin Beta Kernel Process (TwinBKP) model for binary or
#'   binomial response data. TwinBKP is a scalable global-local approximation to
#'   \code{\link{fit_BKP}}. It first selects a global representative subset using
#'   the Twinning algorithm on an augmented representation of the normalized
#'   inputs and empirical response proportions. It then combines a smooth global
#'   kernel contribution from the selected global subset with a compactly
#'   supported local contribution from nearest non-global neighbours.
#'
#' @param X A numeric input matrix of size \eqn{n \times d}, where each row
#'   corresponds to a covariate vector.
#' @param y A numeric vector of observed successes of length \code{n}.
#' @param m A numeric vector of total binomial trials of length \code{n},
#'   corresponding to each \code{y}.
#' @param Xbounds Optional \eqn{d \times 2} matrix specifying the lower and
#'   upper bounds of each input dimension. Used to normalize inputs to
#'   \eqn{[0,1]^d}. If \code{NULL}, inputs are assumed to be pre-normalized, and
#'   default bounds \eqn{[0,1]^d} are applied.
#' @param prior Type of prior: \code{"noninformative"} (default),
#'   \code{"fixed"}, or \code{"adaptive"}.
#' @param r0 Global prior precision (used when \code{prior = "fixed"} or
#'   \code{"adaptive"}).
#' @param p0 Global prior mean (used when \code{prior = "fixed"}). Default is
#'   \code{mean(y/m)}.
#' @param global_kernel Kernel function for the global component:
#'   \code{"gaussian"} (default), \code{"matern52"}, \code{"matern32"}, or
#'   \code{"wendland"}.
#' @param local_kernel Kernel function for the local component. The current
#'   default is \code{"wendland"}, corresponding to the compactly supported
#'   local kernel used by the TwinBKP approximation.
#' @param loss Loss function for global kernel hyperparameter tuning:
#'   \code{"brier"} (default) or \code{"log_loss"}.
#' @param n_multi_start Number of initial points used in multi-start
#'   optimization of the global kernel lengthscale parameters. If \code{NULL},
#'   the default from \code{\link{fit_BKP}} is used on the selected global
#'   subset.
#' @param theta_g Optional. A positive scalar or numeric vector specifying the
#'   global kernel lengthscale parameter(s). If \code{NULL} (default), the global
#'   lengthscale is optimized by fitting a BKP model on the selected global
#'   subset.
#' @param theta_l Optional. A positive scalar specifying the local kernel range.
#'   If \code{NULL} (default), it is set to the empirical covering radius of the
#'   global subset on the normalized input scale.
#' @param isotropic Logical. If \code{TRUE} (default), optimize/use a single
#'   shared global lengthscale across dimensions. If \code{FALSE}, use separate
#'   per-dimension global lengthscales.
#' @param n_threads Number of OpenMP threads used for global hyperparameter
#'   optimization when \code{theta_g = NULL}. Default is \code{1}.
#' @param ess Effective-sample-size calibration for the TwinBKP data
#'   contribution. Use \code{"none"} (default) for the standard update. Use
#'   \code{"shepard"} to rescale the combined global-local data contribution.
#' @param g Target global subset size. If \code{NULL}, the default is
#'   \eqn{\min\{n-1, 50d, \max(\lfloor\sqrt n\rfloor, 10d)\}}. The underlying
#'   C++ Twinning routine follows the \code{twingp} compression-parameter
#'   formulation, so the actual selected size may differ slightly from this
#'   target.
#' @param r Optional Twinning compression parameter. If supplied, it overrides
#'   \code{g}. Larger values produce smaller global subsets.
#' @param l Number of local non-global neighbours used at each training
#'   location. If \code{NULL}, the default is
#'   \eqn{\min\{n-|G|, \max(25, 3d)\}} after the global subset has been selected.
#' @param runs Number of Twinning runs with different starting points.
#' @param u1 Optional integer vector of 1-based starting indices for the
#'   Twinning runs. If \code{NULL}, the first starting point is the observation
#'   farthest from the centroid of the augmented Twinning data, and the remaining
#'   starts are sampled uniformly from the training indices.
#' @param leaf_size Leaf size passed to the \pkg{nanoflann} kd-tree used by the
#'   Twinning routine.
#' @param response_weight Nonnegative scalar weight applied to the empirical
#'   response proportion \code{y/m} when constructing the augmented Twinning
#'   data.
#' @param include_m_in_twin Logical. If \code{TRUE}, a normalized log trial-size
#'   column is included in the augmented Twinning data.
#' @param size_weight Nonnegative scalar weight applied to the normalized log
#'   trial-size column when \code{include_m_in_twin = TRUE}.
#' @param store_kernel Logical. If \code{TRUE}, store dense diagnostic
#'   kernel matrices \code{K}, \code{K_global}, and \code{K_local}. The
#'   default \code{FALSE} avoids \eqn{n \times n} kernel storage and is
#'   required for the stated TwinBKP memory complexity.
#'
#' @return A list of class \code{"TwinBKP"} containing the fitted TwinBKP model,
#'   including:
#' \describe{
#'   \item{\code{theta_opt}}{Alias for the optimized or user-specified global
#'     kernel lengthscale \code{theta_g}.}
#'   \item{\code{theta_g}}{Global kernel lengthscale parameter(s).}
#'   \item{\code{theta_l}}{Local kernel range parameter.}
#'   \item{\code{kernel}}{Alias for \code{global_kernel}, retained for naming
#'     consistency with \code{\link{fit_BKP}}.}
#'   \item{\code{global_kernel}, \code{local_kernel}}{Kernel functions used by
#'     the global and local components.}
#'   \item{\code{isotropic}}{Logical flag indicating whether the global component
#'     uses a shared lengthscale.}
#'   \item{\code{loss}, \code{loss_min}}{Loss function and minimum loss value for
#'     global lengthscale tuning.}
#'   \item{\code{ess}, \code{ess_info}}{ESS calibration method and diagnostics.}
#'   \item{\code{X}, \code{Xnorm}, \code{Xbounds}}{Original inputs, normalized
#'     inputs, and normalization bounds.}
#'   \item{\code{y}, \code{m}}{Observed success counts and binomial trial counts.}
#'   \item{\code{prior}, \code{r0}, \code{p0}}{Prior specification.}
#'   \item{\code{alpha0}, \code{beta0}}{Prior Beta shape parameters at the
#'     training locations.}
#'   \item{\code{alpha_n}, \code{beta_n}}{TwinBKP posterior Beta shape parameters
#'     at the training locations.}
#'   \item{\code{K}, \code{K_global}, \code{K_local}}{Combined, global, and local
#'     kernel-weight matrices used for the training update.}
#'   \item{\code{global_indices}, \code{local_indices}}{Selected global subset and
#'     training-location-specific local subsets, stored as 1-based indices.}
#'   \item{\code{twin_info}, \code{twin_data}}{Diagnostics from the C++ Twinning
#'     routine and the augmented data used for global subset selection.}
#' }
#'
#' @details The global subset is selected using \code{twin_select_global_rcpp()}.
#'   For BKP, the default augmented Twinning data is
#'   \code{cbind(Xnorm, response_weight * y / m)}. Local neighbours are selected
#'   using a kd-tree over non-global training points via \pkg{nanoflann}. By
#'   default, posterior pseudo-counts are aggregated row-wise and dense
#'   \eqn{n \times n} kernel matrices are not stored. Fitting posterior
#'   aggregation is \eqn{O(n(g + l))}. Prediction under \code{ess = "none"}
#'   is \eqn{O(t(\log n + g + l))} for fixed input dimension. Exact
#'   \code{ess = "shepard"} at new prediction points may add \eqn{O(tn)}
#'   arithmetic for exact Shepard interpolation.
#'
#' @seealso \code{\link{fit_BKP}} for the full BKP model,
#'   \code{\link{fit_DKP}} for multinomial responses, \code{\link{predict.BKP}},
#'   \code{\link{plot.BKP}}, \code{\link{simulate.BKP}}, and
#'   \code{\link{summary.BKP}} for downstream workflows on fitted BKP objects.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}. arXiv.
#'
#' Vakayil A, Joseph VR (2022). Data Twinning. \emph{Statistical Analysis and
#'   Data Mining: The ASA Data Science Journal}, 15(5), 598--610.
#'
#' Vakayil A, Joseph VR (2024). A Global-Local Approximation Framework for
#'   Large-Scale Gaussian Process Modeling. \emph{Technometrics}, 66(2),
#'   295--305.
#'
#' @examples
#' set.seed(1)
#' n <- 30
#' X <- matrix(runif(n * 2), ncol = 2)
#' m <- sample(10:30, n, replace = TRUE)
#' p <- plogis(2 * (X[, 1] - 0.5))
#' y <- rbinom(n, size = m, prob = p)
#'
#' model <- fit_TwinBKP(X, y, m, g = 10, runs = 2, theta_g = 0.25)
#' model$theta_g
#' model$theta_l
#'
#' @export
fit_TwinBKP <- function(
    X, y, m, Xbounds = NULL,
    prior = c("noninformative", "fixed", "adaptive"), r0 = 2, p0 = mean(y/m),
    global_kernel = c("gaussian", "matern52", "matern32", "wendland"),
    local_kernel = c("wendland"),
    loss = c("brier", "log_loss"),
    n_multi_start = NULL, theta_g = NULL, theta_l = NULL,
    isotropic = TRUE, n_threads = 1,
    ess = c("none", "shepard"),
    g = NULL, r = NULL, l = NULL,
    runs = 10, u1 = NULL, leaf_size = 8,
    response_weight = 1,
    include_m_in_twin = FALSE,
    size_weight = 0,
    store_kernel = FALSE
) {
  # ---- Argument checking ----
  if (missing(X) || missing(y) || missing(m)) {
    stop("Arguments 'X', 'y', and 'm' must be provided.")
  }
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("'X' must be a numeric matrix or data frame.")
  }
  if (!is.numeric(as.matrix(X))) {
    stop("'X' must contain numeric values only.")
  }
  if (!is.numeric(y)) stop("'y' must be numeric.")
  if (!is.numeric(m)) stop("'m' must be numeric.")

  X <- as.matrix(X)
  y <- as.numeric(y)
  m <- as.numeric(m)

  d <- ncol(X)
  n <- nrow(X)

  if (n < 3L) stop("TwinBKP requires at least three observations.")
  if (length(y) != n) stop("'y' must have the same number of rows as 'X'.")
  if (length(m) != n) stop("'m' must have the same number of rows as 'X'.")
  if (any(y < 0)) stop("'y' must be nonnegative.")
  if (any(m <= 0)) stop("'m' must be strictly positive.")
  if (any(y > m)) stop("Each element of 'y' must be less than or equal to corresponding element of 'm'.")
  if (anyNA(X) || anyNA(y) || anyNA(m)) stop("Missing values are not allowed in 'X', 'y', or 'm'.")
  if (any(!is.finite(X)) || any(!is.finite(y)) || any(!is.finite(m))) {
    stop("'X', 'y', and 'm' must contain only finite values.")
  }

  # ---- prior, kernel, loss ----
  prior <- match.arg(prior)
  global_kernel <- match.arg(global_kernel)
  local_kernel <- match.arg(local_kernel)
  loss <- match.arg(loss)
  ess <- match.arg(ess)

  # ---- Xbounds checks ----
  if (is.null(Xbounds)) {
    xmin <- min(X)
    xmax <- max(X)

    if (xmin < 0 || xmax > 1) {
      warning(
        sprintf(
          paste0(
            "Input X does not appear to be normalized to [0,1]. ",
            "Current range: [%.3f, %.3f].\n",
            "Please normalize X or specify Xbounds explicitly; ",
            "otherwise the model may produce incorrect results."
          ),
          xmin, xmax
        )
      )
    }

    Xbounds <- cbind(rep(0, d), rep(1, d))
  } else {
    if (!is.matrix(Xbounds)) stop("'Xbounds' must be a numeric matrix.")
    if (!is.numeric(Xbounds)) stop("'Xbounds' must contain numeric values.")
    if (!all(dim(Xbounds) == c(d, 2))) {
      stop(paste0("'Xbounds' must be a matrix with dimensions d x 2, where d = ", d, "."))
    }
    if (anyNA(Xbounds) || any(!is.finite(Xbounds))) {
      stop("'Xbounds' must contain only finite values.")
    }
    if (any(Xbounds[,2] <= Xbounds[,1])) {
      stop("Each row of 'Xbounds' must satisfy lower < upper.")
    }
  }

  # ---- prior parameters checks ----
  if (!is.numeric(r0) || length(r0) != 1 || r0 <= 0 || is.na(r0) || !is.finite(r0)) {
    stop("'r0' must be a positive scalar.")
  }

  if (prior == "fixed") {
    if (!is.numeric(p0) || length(p0) != 1 ||
        is.na(p0) || !is.finite(p0) || p0 <= 0 || p0 >= 1) {
      stop("For fixed prior in BKP, 'p0' must be a scalar in (0, 1).")
    }
  }
  p0 <- as.numeric(p0)

  # ---- hyperparameters checks ----
  if (!is.null(n_multi_start)) {
    if (!is.numeric(n_multi_start) || length(n_multi_start) != 1  ||
        is.na(n_multi_start) || !is.finite(n_multi_start) || n_multi_start <= 0) {
      stop("'n_multi_start' must be a positive integer.")
    }
    n_multi_start <- as.integer(n_multi_start)
  }

  if (!is.numeric(n_threads) || length(n_threads) != 1 ||
      is.na(n_threads) || !is.finite(n_threads) || n_threads <= 0) {
    stop("'n_threads' must be a positive integer.")
  }
  n_threads <- as.integer(n_threads)

  if (!is.logical(isotropic) || length(isotropic) != 1) {
    stop("'isotropic' must be a single logical value.")
  }

  if (!is.null(theta_g)) {
    if (!is.numeric(theta_g)) stop("'theta_g' must be numeric.")
    if (isTRUE(isotropic)) {
      if (length(theta_g) != 1) {
        stop("When isotropic=TRUE, 'theta_g' must be a scalar.")
      }
    } else if (!(length(theta_g) == 1 || length(theta_g) == d)) {
      stop(paste0("When isotropic=FALSE, 'theta_g' must be either a scalar or a vector of length ", d, "."))
    }
    if (!isotropic && length(theta_g) == 1) theta_g <- rep(theta_g, d)
    if (anyNA(theta_g) || any(!is.finite(theta_g)) || any(theta_g <= 0)) {
      stop("'theta_g' must be strictly positive.")
    }
    theta_g <- as.numeric(theta_g)
  }

  if (!is.null(theta_l)) {
    if (!is.numeric(theta_l) || length(theta_l) != 1 ||
        is.na(theta_l) || !is.finite(theta_l) || theta_l <= 0) {
      stop("'theta_l' must be a positive scalar.")
    }
    theta_l <- as.numeric(theta_l)
  }

  # ---- Normalize input X to [0,1]^d ----
  Xnorm <- sweep(X, 2, Xbounds[,1], "-")
  Xnorm <- sweep(Xnorm, 2, Xbounds[,2] - Xbounds[,1], "/")

  if (identical(ess, "shepard")) {
    .bkp_check_unique_locations(Xnorm)
  }

  # ---- Twinning controls ----
  if (!is.numeric(runs) || length(runs) != 1 ||
      is.na(runs) || !is.finite(runs) || runs <= 0) {
    stop("'runs' must be a positive integer.")
  }
  runs <- as.integer(runs)

  if (!is.numeric(leaf_size) || length(leaf_size) != 1 ||
      is.na(leaf_size) || !is.finite(leaf_size) || leaf_size <= 0) {
    stop("'leaf_size' must be a positive integer.")
  }
  leaf_size <- as.integer(leaf_size)

  if (!is.numeric(response_weight) || length(response_weight) != 1 ||
      is.na(response_weight) || !is.finite(response_weight) || response_weight < 0) {
    stop("'response_weight' must be a nonnegative scalar.")
  }

  if (!is.logical(include_m_in_twin) || length(include_m_in_twin) != 1) {
    stop("'include_m_in_twin' must be a single logical value.")
  }

  if (!is.logical(store_kernel) || length(store_kernel) != 1) {
    stop("'store_kernel' must be a single logical value.")
  }

  if (!is.numeric(size_weight) || length(size_weight) != 1 ||
      is.na(size_weight) || !is.finite(size_weight) || size_weight < 0) {
    stop("'size_weight' must be a nonnegative scalar.")
  }

  if (is.null(g)) {
    g <- .bkp_twin_default_g(n, d)
  } else {
    if (!is.numeric(g) || length(g) != 1 ||
        is.na(g) || !is.finite(g) || g < 2 || g >= n) {
      stop("'g' must be an integer between 2 and n - 1.")
    }
    g <- as.integer(g)
  }

  if (is.null(r)) {
    r <- ceiling(n / g)
  } else {
    if (!is.numeric(r) || length(r) != 1 ||
        is.na(r) || !is.finite(r) || r < 2 || r >= n) {
      stop("'r' must be an integer between 2 and n - 1.")
    }
    r <- as.integer(r)
  }

  twin_data <- .bkp_build_twin_data_bkp(
    Xnorm = Xnorm,
    y = y,
    m = m,
    response_weight = response_weight,
    include_m_in_twin = include_m_in_twin,
    size_weight = size_weight
  )

  if (is.null(u1)) {
    u1 <- .bkp_twin_starts(twin_data, runs)
  } else {
    if (!is.numeric(u1) || length(u1) != runs || anyNA(u1) || any(!is.finite(u1))) {
      stop("'u1' must be an integer vector of length 'runs'.")
    }
    if (any(u1 < 1) || any(u1 > n)) {
      stop("'u1' must contain valid 1-based row indices.")
    }
    u1 <- as.integer(u1)
  }

  twin_info <- twin_select_global_rcpp(
    twin_data = twin_data,
    Xnorm = Xnorm,
    r = r,
    runs = runs,
    u1 = u1,
    leaf_size = leaf_size
  )

  g_indices <- as.integer(twin_info$g_indices)
  g_actual <- length(g_indices)
  if (g_actual < 2L) {
    stop("The Twinning step selected fewer than two global points; decrease 'r'.")
  }

  if (is.null(theta_l)) {
    theta_l <- as.numeric(twin_info$theta_l)
  }
  if (!is.finite(theta_l) || theta_l <= 0) {
    stop("The local range 'theta_l' must be positive.")
  }

  non_global_n <- n - g_actual
  if (is.null(l)) {
    l <- min(non_global_n, max(25L, 3L * d))
  } else {
    if (!is.numeric(l) || length(l) != 1 ||
        is.na(l) || !is.finite(l) || l < 0) {
      stop("'l' must be a nonnegative integer.")
    }
    l <- as.integer(l)
    if (l > non_global_n) {
      stop("'l' cannot exceed the number of non-global training points.")
    }
  }

  # ---- Global lengthscale tuning on the selected global subset ----
  X_global <- X[g_indices, , drop = FALSE]
  Xnorm_global <- Xnorm[g_indices, , drop = FALSE]
  y_global <- y[g_indices]
  m_global <- m[g_indices]

  if (is.null(theta_g)) {
    global_fit <- fit_BKP(
      X = X_global,
      y = y_global,
      m = m_global,
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

    theta_g <- as.numeric(global_fit$theta_opt)
    loss_min <- as.numeric(global_fit$loss_min)
  } else {
    loss_min <- loss_fun(
      gamma = log10(theta_g),
      Xnorm = Xnorm_global,
      y = y_global,
      m = m_global,
      prior = prior,
      r0 = r0,
      p0 = p0,
      model = "BKP",
      loss = loss,
      kernel = global_kernel,
      isotropic = isotropic,
      ess = ess
    )
  }

  # ---- Compute posterior parameters ----
  local_indices <- twin_local_indices_rcpp(
    Xtrain_norm = Xnorm,
    Xquery_norm = Xnorm,
    g_indices = g_indices,
    l = l,
    leaf_size = leaf_size
  )

  posterior <- .twin_bkp_compute_posterior_fast(
    Xquery_norm = Xnorm,
    Xtrain_norm = Xnorm,
    y = y,
    m = m,
    g_indices = g_indices,
    local_indices = local_indices,
    theta_g = theta_g,
    theta_l = theta_l,
    global_kernel = global_kernel,
    local_kernel = local_kernel,
    isotropic = isotropic,
    prior = prior,
    r0 = r0,
    p0 = p0,
    ess = ess,
    store_kernel = store_kernel
  )

  # ---- Construct and return the fitted model ----
  TwinBKP_model <- list(
    theta_opt = theta_g, theta_g = theta_g, theta_l = theta_l,
    kernel = global_kernel, global_kernel = global_kernel, local_kernel = local_kernel,
    isotropic = isotropic, loss = loss, loss_min = loss_min,
    ess = ess, ess_info = posterior$ess_info,
    X = X, Xnorm = Xnorm, Xbounds = Xbounds, y = matrix(y, ncol = 1), m = matrix(m, ncol = 1),
    prior = prior, r0 = r0, p0 = p0, alpha0 = posterior$alpha0, beta0 = posterior$beta0,
    alpha_n = posterior$alpha_n, beta_n = posterior$beta_n,
    K = posterior$K, K_global = posterior$K_global, K_local = posterior$K_local,
    twin_data = twin_data, twin_info = twin_info,
    global_indices = g_indices, local_indices = local_indices,
    g_target = g, g = g_actual, r = r, l = l, runs = runs, u1 = u1,
    response_weight = response_weight, include_m_in_twin = include_m_in_twin,
    size_weight = size_weight, leaf_size = leaf_size,
    store_kernel = store_kernel,
    complexity = list(
      global_selection = "Twinning global subset selection using kd-tree nearest-neighbour search",
      global_tuning = "O(T_g g^2), where T_g is the number of global-subset loss evaluations",
      training_posterior = "O(n(g + l)) row-wise pseudo-count aggregation",
      prediction = "O(t(log n + g + l)) for ess = 'none' and fixed input dimension",
      memory = "O(n + n l + g) by default, excluding input storage; dense kernels are stored only when store_kernel = TRUE",
      ess_note = "With exact ess = 'shepard', new-point prediction may require O(t n) arithmetic for exact Shepard interpolation."
    )
  )

  class(TwinBKP_model) <- "TwinBKP"
  return(TwinBKP_model)
}


.bkp_twin_default_g <- function(n, d) {
  g <- min(n - 1L, 50L * d, max(floor(sqrt(n)), 10L * d))
  as.integer(max(2L, g))
}


.bkp_build_twin_data_bkp <- function(Xnorm, y, m,
                                     response_weight = 1,
                                     include_m_in_twin = FALSE,
                                     size_weight = 0) {
  p_hat <- as.numeric(y) / as.numeric(m)
  twin_data <- cbind(Xnorm, response_weight * p_hat)

  if (isTRUE(include_m_in_twin)) {
    log_m <- log1p(as.numeric(m))
    rng <- range(log_m)

    if (diff(rng) > 0) {
      m_scaled <- (log_m - rng[1L]) / diff(rng)
    } else {
      m_scaled <- rep(0, length(log_m))
    }

    twin_data <- cbind(twin_data, size_weight * m_scaled)
  }

  storage.mode(twin_data) <- "double"
  twin_data
}


.bkp_twin_starts <- function(twin_data, runs) {
  n <- nrow(twin_data)
  center <- colMeans(twin_data)
  d2 <- rowSums(sweep(twin_data, 2, center, "-")^2)

  if (runs <= 1L) {
    return(as.integer(which.max(d2)))
  }

  as.integer(c(which.max(d2), sample.int(n, runs - 1L, replace = TRUE)))
}


.bkp_twin_local_indices <- function(Xnorm, g_indices, l) {
  n <- nrow(Xnorm)

  if (l == 0L) {
    return(matrix(integer(0), nrow = n, ncol = 0L))
  }

  non_global <- setdiff(seq_len(n), as.integer(g_indices))
  out <- matrix(NA_integer_, nrow = n, ncol = l)

  for (i in seq_len(n)) {
    dif <- sweep(Xnorm[non_global, , drop = FALSE], 2, Xnorm[i, ], "-")
    d2 <- rowSums(dif^2)
    out[i, ] <- non_global[order(d2)[seq_len(l)]]
  }

  out
}


.twin_bkp_compute_posterior <- function(Xnorm, y, m, g_indices, local_indices,
                                        theta_g, theta_l,
                                        global_kernel, local_kernel, isotropic,
                                        prior, r0, p0, ess = "none") {
  n <- nrow(Xnorm)
  y <- as.numeric(y)
  m <- as.numeric(m)
  g_indices <- as.integer(g_indices)

  K_global <- matrix(0, nrow = n, ncol = n)
  K_global[, g_indices] <- kernel_matrix(
    X = Xnorm,
    Xprime = Xnorm[g_indices, , drop = FALSE],
    theta = theta_g,
    kernel = global_kernel,
    isotropic = isotropic
  )

  K_local <- matrix(0, nrow = n, ncol = n)
  if (ncol(local_indices) > 0L) {
    for (i in seq_len(n)) {
      idx <- local_indices[i, ]
      K_local[i, idx] <- as.numeric(kernel_matrix(
        X = Xnorm[i, , drop = FALSE],
        Xprime = Xnorm[idx, , drop = FALSE],
        theta = theta_l,
        kernel = local_kernel,
        isotropic = TRUE
      ))
    }
  }

  K <- K_global + K_local

  prior_par <- get_prior(
    prior = prior,
    model = "BKP",
    r0 = r0,
    p0 = p0,
    y = y,
    m = m,
    K = K
  )
  alpha0 <- prior_par$alpha0
  beta0 <- prior_par$beta0

  data_success <- as.vector(K %*% y)
  data_failure <- as.vector(K %*% (m - y))

  if (identical(ess, "shepard")) {
    ess_info <- .bkp_ess_calibration(
      Xquery_norm = Xnorm,
      Xtrain_norm = Xnorm,
      m = m,
      K = K
    )
    data_success <- ess_info$scale * data_success
    data_failure <- ess_info$scale * data_failure
  } else {
    ess_info <- .bkp_ess_none_info(K, m)
  }

  list(
    K = K,
    K_global = K_global,
    K_local = K_local,
    alpha0 = alpha0,
    beta0 = beta0,
    alpha_n = alpha0 + data_success,
    beta_n = beta0 + data_failure,
    ess_info = ess_info
  )
}

.twin_bkp_compute_posterior_fast <- function(
    Xquery_norm, Xtrain_norm, y, m, g_indices, local_indices,
    theta_g, theta_l, global_kernel, local_kernel, isotropic,
    prior, r0, p0, ess = "none", store_kernel = FALSE
) {
  m_shepard <- NULL

  if (identical(ess, "shepard")) {
    m_shepard <- .bkp_shepard_m(
      Xquery_norm = Xquery_norm,
      Xtrain_norm = Xtrain_norm,
      m = as.numeric(m),
      power = 2
    )
  }

  twin_bkp_posterior_rcpp(
    Xquery_norm = Xquery_norm,
    Xtrain_norm = Xtrain_norm,
    y = as.numeric(y),
    m = as.numeric(m),
    g_indices = as.integer(g_indices),
    local_indices = local_indices,
    theta_g = as.numeric(theta_g),
    theta_l = as.numeric(theta_l),
    global_kernel = global_kernel,
    local_kernel = local_kernel,
    isotropic = isTRUE(isotropic),
    prior = prior,
    r0 = r0,
    p0 = p0,
    ess = ess,
    m_shepard = m_shepard,
    store_kernel = store_kernel
  )
}
