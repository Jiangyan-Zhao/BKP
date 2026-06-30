#' @title Leave-One-Out Loss for BKP and DKP Models
#'
#' @description Computes the leave-one-out cross-validation loss used for
#'   kernel hyperparameter tuning in BKP and DKP models. The input
#'   \code{gamma} represents log10-transformed kernel lengthscale parameters,
#'   with \code{theta = 10^gamma}. The function supports Brier score and
#'   log-loss under noninformative, fixed, and adaptive prior specifications.
#'
#' @inheritParams get_prior
#' @inheritParams fit_BKP
#' @param gamma A numeric vector of log10-transformed kernel hyperparameters.
#' @param Xnorm A numeric matrix of normalized input features (\code{[0,1]^d}).
#' @param ess Effective-sample-size calibration for leave-one-out data
#'   contributions. Use \code{"none"} (default) for the standard LOOCV loss.
#'   Use \code{"shepard"} to apply the same data-precision rescaling as in
#'   model fitting, with the Shepard target computed in leave-one-out form.
#'
#' @return A single numeric value giving the average leave-one-out loss to be
#'   minimized. For BKP, the Brier score and log-loss are averaged over
#'   observations. For DKP, class-wise losses are summed within each observation
#'   and then averaged over observations.
#'
#' @details The kernel lengthscale parameter is represented on the log10 scale:
#'   \code{theta = 10^gamma}. The function first computes the kernel matrix
#'   using \code{\link{kernel_matrix}} and then sets its diagonal entries to
#'   zero, so that each fitted value is computed with the corresponding
#'   observation excluded. This implements leave-one-out cross-validation
#'   without repeatedly refitting the model.
#'
#'   For \code{model = "BKP"}, the loss compares the leave-one-out posterior
#'   mean success probability with the empirical proportion \code{y / m}. For
#'   \code{model = "DKP"}, the loss compares the leave-one-out posterior mean
#'   class-probability vector with the empirical class-proportion vector
#'   \code{Y / rowSums(Y)}.
#'
#'   If \code{ess = "shepard"}, the same data-precision rescaling used in model
#'   fitting is applied in leave-one-out form. This option requires unique
#'   training input locations on the normalized input scale.
#'
#' @seealso \code{\link{fit_BKP}} for fitting BKP models, \code{\link{fit_DKP}}
#'   for fitting DKP models, \code{\link{get_prior}} for constructing prior
#'   parameters, \code{\link{kernel_matrix}} for computing kernel matrices.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv. <doi:10.48550/arXiv.2508.10447>.
#'
#' @examples
#' # -------------------------- BKP ---------------------------
#' set.seed(123)
#'
#' n <- 10
#' Xnorm <- matrix(runif(2 * n), ncol = 2)
#' m <- rep(10, n)
#' pi_true <- runif(n)
#' y <- rbinom(n, size = m, prob = pi_true)
#' loss_bkp <- loss_fun(gamma = log10(0.2), Xnorm = Xnorm, y = y, m = m)
#'
#' # -------------------------- DKP ---------------------------
#' n <- 10
#' q <- 3
#' Xnorm <- matrix(runif(2 * n), ncol = 2)
#' m <- rep(10, n)
#' p0 <- rep(1 / q, q)
#' Y <- matrix(rmultinom(n, size = 10, prob = rep(1/q, q)), nrow = n, byrow = TRUE)
#' loss_dkp <- loss_fun(gamma = log10(0.2), Xnorm = Xnorm, Y = Y, model = "DKP")
#'
#' @export
loss_fun <- function(
    gamma, Xnorm,
    y = NULL, m = NULL, Y = NULL,
    model = c("BKP", "DKP"),
    prior = c("noninformative", "fixed", "adaptive"), r0 = 2, p0 = NULL,
    loss = c("brier", "log_loss"),
    kernel = c("gaussian", "matern52", "matern32", "wendland"),
    isotropic = TRUE,
    ess = c("none", "shepard"))
{
  # ---- Argument checking ----
  if (!is.numeric(gamma) || length(gamma) < 1L ||
      anyNA(gamma) || any(!is.finite(gamma))) {
    stop("'gamma' must be a finite numeric vector.")
  }

  if (!is.matrix(Xnorm) || !is.numeric(Xnorm) ||
      anyNA(Xnorm) || any(!is.finite(Xnorm))) {
    stop("'Xnorm' must be a numeric matrix with finite values and no NA values.")
  }
  if (nrow(Xnorm) < 1L || ncol(Xnorm) < 1L) {
    stop("'Xnorm' must have at least one row and one column.")
  }

  model <- match.arg(model)
  prior <- match.arg(prior)
  loss <- match.arg(loss)
  kernel <- match.arg(kernel)
  ess <- match.arg(ess)

  if (!is.numeric(r0) || length(r0) != 1L ||
      is.na(r0) || !is.finite(r0) || r0 <= 0) {
    stop("'r0' must be a positive finite scalar.")
  }

  if (!is.null(p0)) {
    if (!is.numeric(p0) || anyNA(p0) || any(!is.finite(p0)) || any(p0 < 0)) {
      stop("'p0' must contain nonnegative finite numeric values.")
    }
  }

  if (!is.logical(isotropic) || length(isotropic) != 1L) {
    stop("'isotropic' must be a single logical value.")
  }

  if (model == "BKP") {
    if (is.null(y) || is.null(m)) {
      stop("'y' and 'm' must be provided for BKP model.")
    }
    if (!is.numeric(y) || !is.numeric(m) ||
        anyNA(y) || anyNA(m) ||
        any(!is.finite(y)) || any(!is.finite(m))) {
      stop("'y' and 'm' must be numeric vectors with finite values and no NA values.")
    }

    y <- as.numeric(y)
    m <- as.numeric(m)

    if (length(y) != nrow(Xnorm) || length(m) != nrow(Xnorm)) {
      stop("'y' and 'm' must have the same length as the number of rows in 'Xnorm'.")
    }
    if (any(y < 0) || any(m <= 0) || any(y > m)) {
      stop("'y' must be in [0, m] and 'm' must be positive.")
    }
  } else {
    if (is.null(Y)) {
      stop("'Y' must be provided for DKP model.")
    }
    if (!is.matrix(Y) || !is.numeric(Y) ||
        anyNA(Y) || any(!is.finite(Y))) {
      stop("'Y' must be a numeric matrix with finite values and no NA values.")
    }
    if (nrow(Y) != nrow(Xnorm)) {
      stop("Number of rows in 'Y' must match the number of rows in 'Xnorm'.")
    }
    if (ncol(Y) < 2L) {
      stop("'Y' must have at least two columns.")
    }
    if (any(Y < 0)) {
      stop("'Y' must contain nonnegative counts or frequencies.")
    }
    if (any(rowSums(Y) <= 0)) {
      stop("Each row of 'Y' must have a positive row sum.")
    }
  }

  if (identical(ess, "shepard")) {
    bkp_check_unique_locations(Xnorm)
    if (nrow(Xnorm) < 2L) {
      stop("ESS calibration with ess = 'shepard' requires at least two input locations for leave-one-out calibration.")
    }
  }

  # Convert gamma to kernel hyperparameters (theta = 10^gamma)
  theta <- 10^gamma

  # Compute kernel matrix using specified kernel and theta
  K <- kernel_matrix(Xnorm, theta = theta, kernel = kernel, isotropic = isotropic)

  diag(K) <- 0  # Leave-One-Out Cross-Validation (LOOCV)

  if (model == "BKP") {
    # Get prior parameters
    prior_par <- get_prior(prior = prior, model = model, r0 = r0, p0 = p0, y = y, m = m, K = K)

    data_scale <- NULL
    if (identical(ess, "shepard")) {
      m_kernel <- as.vector(K %*% as.numeric(m))
      rho <- apply(K, 1L, max)
      m_target <- rho * bkp_shepard_m_loo(Xnorm, m, power = 2)
      data_scale <- numeric(length(m_kernel))
      positive_kernel_mass <- m_kernel > 0
      data_scale[positive_kernel_mass] <- m_target[positive_kernel_mass] / m_kernel[positive_kernel_mass]
    }

    result <- loss_fun_rcpp(
      model = model,
      loss = loss,
      K = K,
      y = as.numeric(y),
      m = as.numeric(m),
      Y = NULL,
      alpha0 = as.numeric(prior_par$alpha0),
      beta0 = as.numeric(prior_par$beta0),
      alpha0_mat = NULL,
      data_scale = data_scale
    )
  } else {
    # Get prior parameters for DKP
    alpha0 <- get_prior(prior = prior, model = model, r0 = r0, p0 = p0, Y = Y, K = K)

    data_scale <- NULL
    if (ess == "shepard") {
      m_dkp <- rowSums(Y)
      m_kernel <- as.vector(K %*% as.numeric(m_dkp))
      rho <- apply(K, 1L, max)
      m_target <- rho * bkp_shepard_m_loo(Xnorm, m_dkp, power = 2)
      data_scale <- numeric(length(m_kernel))
      positive_kernel_mass <- m_kernel > 0
      data_scale[positive_kernel_mass] <- m_target[positive_kernel_mass] / m_kernel[positive_kernel_mass]
    }

    result <- loss_fun_rcpp(
      model = model,
      loss = loss,
      K = K,
      y = NULL,
      m = NULL,
      Y = as.matrix(Y),
      alpha0 = NULL,
      beta0 = NULL,
      alpha0_mat = as.matrix(alpha0),
      data_scale = data_scale
    )
  }

  return(result)
}
