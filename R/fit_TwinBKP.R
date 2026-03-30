#' @name fit_TwinBKP
#'
#' @title Fit a Twin Beta Kernel Process (TwinBKP) Model
#'
#' @description
#' Fits a TwinBKP model for binary or binomial response data.
#' The workflow follows \code{fit_BKP()}, but kernel hyperparameter tuning
#' is accelerated by first selecting \code{g} representative support points
#' from the full dataset using the \code{Twining} package, then optimizing
#' the kernel lengthscales on these \code{g} points only.
#' The final posterior update still uses all \code{n} observations.
#'
#' @inheritParams fit_BKP
#' @param g_nums Positive integer. Number of global support points selected by
#'   the Twining algorithm for hyperparameter tuning. If \code{NULL} (default),
#'   it is set adaptively as
#'   \eqn{\min\{50d, \max(10d, \sqrt{n})\}},
#'   where \eqn{d} is input dimensionality and \eqn{n} is sample size.
#'
#' @return A list of class \code{"TwinBKP"} with the same structure as
#'   \code{"BKP"}, plus:
#' \describe{
#'   \item{\code{g_nums}}{Number of support points used for tuning.}
#'   \item{\code{tune_idx}}{Row indices of the selected support points.}
#' }
#'
#' @seealso \code{\link{fit_BKP}}
#'
#' @examples
#' set.seed(123)
#' true_pi_fun <- function(x) {
#'   (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' }
#' n <- 50
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#' model <- fit_TwinBKP(X, y, m, Xbounds = Xbounds, g_nums = 20)
#' print(model)
#'
#' @export
fit_TwinBKP <- function(
    X, y, m, Xbounds = NULL,
    prior = c("noninformative", "fixed", "adaptive"), r0 = 2, p0 = mean(y / m),
    kernel = c("gaussian", "matern52", "matern32"),
    loss = c("brier", "log_loss"),
    n_multi_start = NULL, theta = NULL,
    isotropic = TRUE,
    g_nums = NULL
) { 

  if (missing(X) || missing(y) || missing(m)) {
    stop("Arguments 'X', 'y', and 'm' must be provided.")
  }
  if (!is.matrix(X) && !is.data.frame(X)) stop("'X' must be a numeric matrix or data frame.")
  if (!is.numeric(as.matrix(X))) stop("'X' must contain numeric values only.")
  if (!is.numeric(y)) stop("'y' must be numeric.")
  if (!is.numeric(m)) stop("'m' must be numeric.") 
  
  X <- as.matrix(X)
  y <- matrix(y, ncol = 1)
  m <- matrix(m, ncol = 1)

  d <- ncol(X)
  n <- nrow(X)

  if (nrow(y) != n) stop("'y' must have the same number of rows as 'X'.")
  if (nrow(m) != n) stop("'m' must have the same number of rows as 'X'.")

  if (any(y < 0))  stop("'y' must be nonnegative.")
  if (any(m <= 0)) stop("'m' must be strictly positive.")
  if (any(y > m))  stop("Each element of 'y' must be <= corresponding element of 'm'.")
  if (anyNA(X) || anyNA(y) || anyNA(m)) {
    stop("Missing values are not allowed in 'X', 'y', or 'm'.")
  }

  prior  <- match.arg(prior)
  kernel <- match.arg(kernel)
  loss   <- match.arg(loss)

  if (is.null(g_nums)) {
    g_default <- min(50 * d, max(10 * d, sqrt(n)))
    g_nums <- as.integer(round(g_default))
  }
  if (!is.numeric(g_nums) || length(g_nums) != 1 || g_nums <= 0) {
    stop("'g_nums' must be a positive integer.")
  }
  g_nums <- as.integer(g_nums)
  g_eff <- min(g_nums, n)

  if (is.null(Xbounds)) {
    xmin <- min(X); xmax <- max(X)
    if (xmin < 0 || xmax > 1) {
      warning(sprintf(
        paste0(
          "Input X does not appear to be normalized to [0,1]. ",
          "Current range: [%.3f, %.3f]. ",
          "Please specify Xbounds explicitly."
        ), xmin, xmax
      ))
    }
    Xbounds <- cbind(rep(0, d), rep(1, d))
  } else {
    if (!is.matrix(Xbounds))   stop("'Xbounds' must be a numeric matrix.")
    if (!is.numeric(Xbounds))  stop("'Xbounds' must contain numeric values.")
    if (!all(dim(Xbounds) == c(d, 2))) {
      stop(paste0("'Xbounds' must be a d x 2 matrix, where d = ", d, "."))
    }
    if (any(Xbounds[, 2] <= Xbounds[, 1])) {
      stop("Each row of 'Xbounds' must satisfy lower < upper.")
    }
  }

  Xnorm <- sweep(X, 2, Xbounds[, 1], "-")
  Xnorm <- sweep(Xnorm, 2, Xbounds[, 2] - Xbounds[, 1], "/")


  if (!is.numeric(r0) || length(r0) != 1 || r0 <= 0) stop("'r0' must be a positive scalar.")
  if (!is.numeric(p0) || length(p0) != 1 || p0 <= 0 || p0 >= 1) stop("'p0' must be in (0,1).")

  if (!is.null(n_multi_start)) {
    if (!is.numeric(n_multi_start) || length(n_multi_start) != 1 || n_multi_start <= 0) {
      stop("'n_multi_start' must be a positive integer.")
    }
  }

  if (!is.logical(isotropic) || length(isotropic) != 1) {
    stop("'isotropic' must be a single logical value.")
  }

  if (!is.null(theta)) {
    if (!is.numeric(theta)) stop("'theta' must be numeric.")
    if (isotropic && length(theta) != 1) {
      stop("When isotropic=TRUE, 'theta' must be a scalar.")
    }
    if (!isotropic && !(length(theta) == 1 || length(theta) == d)) {
      stop(paste0("When isotropic=FALSE, 'theta' must be scalar or length ", d, "."))
    }
    if (!isotropic && length(theta) == 1) theta <- rep(theta, d)
    if (any(theta <= 0)) stop("'theta' must be strictly positive.")
  }

  tw <- get_twin_indices(Xnorm, g = g_eff, v = 2L * g_eff, runs = 10L, seed = 123L)
  global_idx <- as.integer(tw$gIndices)

  X_global    <- X[global_idx, , drop = FALSE]
  Xnorm_global <- Xnorm[global_idx, , drop = FALSE]
  y_global    <- y[global_idx, , drop = FALSE]
  m_global    <- m[global_idx, , drop = FALSE]


  if (is.null(theta)) {
    n_theta <- if (isotropic) 1L else d

    n_grid_cpp <- as.integer(max(10L, 10L * d))

    n_starts_cpp <- if (is.null(n_multi_start)) {
      1L
    } else {
      as.integer(max(1L, n_multi_start))
    }

    max_iter_cpp <- 100L 
    g_lower <- (log10(d) - log10(500)) / 2
    g_upper <- (log10(d) + 2) / 2
    opt_cpp <- optimize_bkp_theta_rcpp(
      Xnorm = Xnorm_global,
      y = as.numeric(y_global),
      m = as.numeric(m_global),
      prior = prior,
      r0 = r0,
      p0 = p0,
      loss = loss,
      kernel = kernel,
      isotropic = isotropic,
      n_grid = n_grid_cpp,
      n_starts = n_starts_cpp,
      max_iter = max_iter_cpp,
      g_lower = g_lower,
      g_upper = g_upper
    )

    gamma_opt    <- as.numeric(opt_cpp$gamma_opt)
    theta_global <- as.numeric(opt_cpp$theta_opt)
    loss_global  <- as.numeric(opt_cpp$loss_min)
  } else {
    theta_global <- theta
    loss_global  <- loss_fun(
      gamma  = log10(theta_global),
      Xnorm  = Xnorm_global,
      y      = y_global,
      m      = m_global,
      prior  = prior, r0 = r0, p0 = p0,
      model  = "BKP", loss = loss,
      kernel = kernel, isotropic = isotropic
    )
  }


  TwinBKP_model <- list(

    theta_global = theta_global,
    loss_global  = loss_global, 

    global_idx   = global_idx,
    X_global     = X_global,
    Xnorm_global = Xnorm_global,
    y_global     = y_global,
    m_global     = m_global,

    X      = X,
    Xnorm  = Xnorm,
    Xbounds = Xbounds,
    y      = y,
    m      = m,

    kernel    = kernel,
    isotropic = isotropic,
    prior     = prior,
    r0        = r0,
    p0        = p0,
    loss      = loss,
    g_nums    = g_eff,

    theta_opt = theta_global,
    loss_min = loss_global,
    tune_idx = global_idx
  )

  class(TwinBKP_model) <- "TwinBKP"
  return(TwinBKP_model)
}