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
    isotropic = TRUE)
{
  # ---- Argument checking ----
  if (!is.numeric(gamma)) stop("'gamma' must be a numeric vector.")
  if (!is.matrix(Xnorm) || anyNA(Xnorm)) stop("'Xnorm' must be a numeric matrix with no NA.")

  model <- match.arg(model)
  prior <- match.arg(prior)
  loss <- match.arg(loss)
  kernel <- match.arg(kernel)

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
    alpha0 <- prior_par$alpha0
    beta0 <- prior_par$beta0

    # Call C++ function
    if (loss == "brier") {
      result <- loss_fun_brier_bkp_rcpp(K, as.numeric(y), as.numeric(m), as.numeric(alpha0), as.numeric(beta0))
    } else {
      result <- loss_fun_logloss_bkp_rcpp(K, as.numeric(y), as.numeric(m), as.numeric(alpha0), as.numeric(beta0))
    }
  } else {
    # Get prior parameters for DKP
    alpha0 <- get_prior(prior = prior, model = model, r0 = r0, p0 = p0, Y = Y, K = K)

    # Call C++ function
    if (loss == "brier") {
      result <- loss_fun_brier_dkp_rcpp(K, as.matrix(Y), as.matrix(alpha0))
    } else {
      result <- loss_fun_logloss_dkp_rcpp(K, as.matrix(Y), as.matrix(alpha0))
    }
  }

  result
}