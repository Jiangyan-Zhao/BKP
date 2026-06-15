#' @title Loss Function for BKP and DKP Models
#'
#' @description Computes the loss for fitting BKP (binary) or DKP (multi-class)
#'   models. Supports Brier score (mean squared error) and log-loss
#'   (cross-entropy) under different prior specifications.
#'
#' @inheritParams get_prior
#' @inheritParams fit_BKP
#' @param gamma A numeric vector of log10-transformed kernel hyperparameters.
#' @param Xnorm A numeric matrix of normalized input features (\code{[0,1]^d}).
#' @param ess Effective-sample-size calibration for leave-one-out data
#'   contributions. Use \code{"none"} (default) for the standard loss or
#'   \code{"shepard"} for leave-one-out Shepard ESS calibration.
#'
#' @return A single numeric value representing the total loss (to be minimized).
#'   The value corresponds to either the Brier score (squared error) or the
#'   log-loss (cross-entropy).
#'
#' @seealso \code{\link{fit_BKP}} for fitting BKP models, \code{\link{fit_DKP}}
#'   for fitting DKP models, \code{\link{get_prior}} for constructing prior
#'   parameters, \code{\link{kernel_matrix}} for computing kernel matrices.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv.
#'   https://doi.org/10.48550/arXiv.2508.10447.
#'
#' @examples
#' # -------------------------- BKP ---------------------------
#' set.seed(123)
#' n <- 10
#' Xnorm <- matrix(runif(2 * n), ncol = 2)
#' m <- rep(10, n)
#' y <- rbinom(n, size = m, prob = runif(n))
#' loss_fun(gamma = 0, Xnorm = Xnorm, y = y, m = m, model = "BKP")
#'
#' # -------------------------- DKP ---------------------------
#' set.seed(123)
#' n <- 10
#' q <- 3
#' Xnorm <- matrix(runif(2 * n), ncol = 2)
#' Y <- matrix(rmultinom(n, size = 10, prob = rep(1/q, q)), nrow = n, byrow = TRUE)
#' loss_fun(gamma = 0, Xnorm = Xnorm, Y = Y, model = "DKP")
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
  if (!is.numeric(gamma)) stop("'gamma' must be a numeric vector.")
  if (!is.matrix(Xnorm) || anyNA(Xnorm)) stop("'Xnorm' must be a numeric matrix with no NA.")

  model <- match.arg(model)
  prior <- match.arg(prior)
  loss <- match.arg(loss)
  kernel <- match.arg(kernel)
  ess <- match.arg(ess)

  if (identical(ess, "shepard")) {
    .bkp_check_unique_locations(Xnorm)
    if (nrow(Xnorm) < 2L) {
      stop("ESS calibration with ess = 'shepard' requires at least two input locations for leave-one-out calibration.")
    }
  }

  if (model == "BKP") {
    if (is.null(y) || is.null(m)) stop("'y' and 'm' must be provided for BKP model.")
    if (!is.numeric(y) || !is.numeric(m)) stop("'y' and 'm' must be numeric vectors.")
    if (any(y < 0) || any(m <= 0) || any(y > m)) stop("'y' must be in [0,m] and 'm' > 0.")
    if (length(y) != nrow(Xnorm) || length(m) != nrow(Xnorm)) {
      stop("'y' and 'm' must have the same length as number of rows in 'Xnorm'.")
    }
  } else {
    if (is.null(Y)) stop("'Y' must be provided for DKP model.")
    if (!is.matrix(Y) || anyNA(Y) || any(Y < 0)) stop("'Y' must be a numeric matrix with no NA and nonnegative entries.")
    if (nrow(Y) != nrow(Xnorm)) stop("Number of rows in 'Y' must match number of rows in 'Xnorm'.")
    if (any(rowSums(Y) <= 0)) stop("Each row of 'Y' must have positive total count.")
  }

  if (!is.numeric(r0) || length(r0) != 1 || r0 <= 0) stop("'r0' must be a positive scalar.")
  if (!is.null(p0) && (!is.numeric(p0) || any(p0 < 0))) stop("'p0' must be numeric and nonnegative.")
  if (!is.logical(isotropic) || length(isotropic) != 1) stop("'isotropic' must be a single logical value.")

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
      m_target <- rho * .bkp_shepard_m_loo(Xnorm, m, power = 2)
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
    if (identical(ess, "shepard")) {
      m_dkp <- rowSums(Y)
      m_kernel <- as.vector(K %*% as.numeric(m_dkp))
      rho <- apply(K, 1L, max)
      m_target <- rho * .bkp_shepard_m_loo(Xnorm, m_dkp, power = 2)
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
