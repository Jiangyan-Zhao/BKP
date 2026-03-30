#' @name fit_TwinDKP
#'
#' @title Fit a Twin Dirichlet Kernel Process (TwinDKP) Model
#'
#' @description
#' Fits a TwinDKP model for categorical or multinomial response data.
#' Similar to \code{fit_TwinBKP()}, the method first selects \code{g_nums}
#' representative support points using Twining, then optimizes kernel
#' hyperparameters on these support points only.
#'
#' @inheritParams fit_DKP
#' @param g_nums Positive integer. Number of global support points selected by
#'   the Twining algorithm for hyperparameter tuning. If \code{NULL} (default),
#'   it is set adaptively as
#'   \eqn{\min\{50d, \max(10d, \sqrt{n})\}},
#'   where \eqn{d} is input dimensionality and \eqn{n} is sample size.
#'
#' @return A list of class \code{"TwinDKP"} with the same structure as
#'   \code{"DKP"}, plus:
#' \describe{
#'   \item{\code{g_nums}}{Number of support points used for tuning.}
#'   \item{\code{tune_idx}}{Row indices of the selected support points.}
#' }
#'
#' @seealso \code{\link{fit_DKP}}, \code{\link{predict.DKP}}
#' @export
fit_TwinDKP <- function(
    X, Y, Xbounds = NULL,
    prior = c("noninformative", "fixed", "adaptive"), r0 = 2, p0 = colMeans(Y / rowSums(Y)),
    kernel = c("gaussian", "matern52", "matern32"),
    loss = c("brier", "log_loss"),
    n_multi_start = NULL, theta = NULL,
    isotropic = TRUE,
    g_nums = NULL
) {
  if (missing(X) || missing(Y)) {
    stop("Arguments 'X' and 'Y' must be provided.")
  }
  if (!is.matrix(X) && !is.data.frame(X)) stop("'X' must be a numeric matrix or data frame.")
  if (!is.numeric(as.matrix(X))) stop("'X' must contain numeric values only.")
  if (!is.matrix(Y) && !is.data.frame(Y)) stop("'Y' must be a numeric matrix or data frame.")
  if (!is.numeric(as.matrix(Y))) stop("'Y' must contain numeric values only.")

  X <- as.matrix(X)
  Y <- as.matrix(Y)

  d <- ncol(X)
  n <- nrow(X)
  q <- ncol(Y)

  if (nrow(Y) != n) stop("Number of rows in 'Y' must match number of rows in 'X'.")
  if (any(Y < 0)) stop("'Y' must be nonnegative counts or frequencies.")
  if (anyNA(X) || anyNA(Y)) stop("Missing values are not allowed in 'X' or 'Y'.")
  if (q < 2) stop("'Y' must have at least two columns (multinomial outcomes).")

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
    if (!is.matrix(Xbounds)) stop("'Xbounds' must be a numeric matrix.")
    if (!is.numeric(Xbounds)) stop("'Xbounds' must contain numeric values.")
    if (!all(dim(Xbounds) == c(d, 2))) {
      stop(paste0("'Xbounds' must be a d x 2 matrix, where d = ", d, "."))
    }
    if (any(Xbounds[, 2] <= Xbounds[, 1])) {
      stop("Each row of 'Xbounds' must satisfy lower < upper.")
    }
  }

  Xnorm <- sweep(X, 2, Xbounds[, 1], "-")
  Xnorm <- sweep(Xnorm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

  if (!is.numeric(r0) || length(r0) != 1 || r0 <= 0) {
    stop("'r0' must be a positive scalar.")
  }
  if (!is.numeric(p0) || any(p0 < 0) || abs(sum(p0) - 1) > 1e-10) {
    stop("'p0' must be numeric, nonnegative, and sum to 1.")
  }

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
    if (isotropic && length(theta) != 1) stop("When isotropic=TRUE, 'theta' must be a scalar.")
    if (!isotropic && !(length(theta) == 1 || length(theta) == d)) {
      stop(paste0("When isotropic=FALSE, 'theta' must be scalar or length ", d, "."))
    }
    if (!isotropic && length(theta) == 1) theta <- rep(theta, d)
    if (any(theta <= 0)) stop("'theta' must be strictly positive.")
  }

  tw <- get_twin_indices(Xnorm, g = g_eff, v = 2L * g_eff, runs = 10L, seed = 123L)
  global_idx <- as.integer(tw$gIndices)

  X_global <- X[global_idx, , drop = FALSE]
  Xnorm_global <- Xnorm[global_idx, , drop = FALSE]
  Y_global <- Y[global_idx, , drop = FALSE]

  if (is.null(theta)) {
    n_theta <- if (isotropic) 1 else d
    gamma_bounds <- matrix(c((log10(d)-log10(500))/2,
                             (log10(d)+2)/2),
                           ncol = 2, nrow = n_theta, byrow = TRUE)
    if (is.null(n_multi_start)) n_multi_start <- 10 * n_theta
    init_gamma <- tgp::lhs(n_multi_start, gamma_bounds)

    opt_res <- optimx::multistart(
      parmat = init_gamma,
      fn     = loss_fun,
      method = "L-BFGS-B",
      lower  = rep(-3, n_theta),
      upper  = rep(3, n_theta),
      prior = prior, r0 = r0, p0 = p0,
      Xnorm = Xnorm_global, Y = Y_global,
      model = "DKP", loss = loss, kernel = kernel, isotropic = isotropic,
      control = list(trace = 0)
    )

    best_index <- which.min(opt_res$value)
    gamma_opt <- as.numeric(opt_res[best_index, 1:n_theta])
    theta_global <- 10^gamma_opt
    loss_global <- opt_res$value[best_index]
  } else {
    theta_global <- theta
    loss_global <- loss_fun(
      gamma = log10(theta_global),
      Xnorm = Xnorm_global, Y = Y_global,
      prior = prior, r0 = r0, p0 = p0,
      model = "DKP", loss = loss, kernel = kernel,
      isotropic = isotropic
    )
  }

  TwinDKP_model <- list(
    theta_global = theta_global,
    loss_global  = loss_global,

    global_idx   = global_idx,
    X_global     = X_global,
    Xnorm_global = Xnorm_global,
    Y_global     = Y_global,

    X = X,
    Xnorm = Xnorm,
    Xbounds = Xbounds,
    Y = Y,

    kernel = kernel,
    isotropic = isotropic,
    prior = prior,
    r0 = r0,
    p0 = p0,
    loss = loss,
    g_nums = g_eff,

    theta_opt = theta_global,
    loss_min = loss_global,
    tune_idx = global_idx
  )

  class(TwinDKP_model) <- "TwinDKP"
  TwinDKP_model
}
