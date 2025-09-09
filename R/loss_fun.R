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
#' loss_fun(gamma = rep(0, 2), Xnorm = Xnorm, y = y, m = m, model = "BKP")
#'
#' # -------------------------- DKP ---------------------------
#' set.seed(123)
#' n <- 10
#' q <- 3
#' Xnorm <- matrix(runif(2 * n), ncol = 2)
#' Y <- matrix(rmultinom(n, size = 10, prob = rep(1/q, q)), nrow = n, byrow = TRUE)
#' loss_fun(gamma = rep(0, 2), Xnorm = Xnorm, Y = Y, model = "DKP")
#'
#' @export

loss_fun <- function(
    gamma, Xnorm,
    y = NULL, m = NULL, Y = NULL,
    model = c("BKP", "DKP"),
    prior = c("noninformative", "fixed", "adaptive"), r0 = 2, p0 = NULL,
    loss = c("brier", "log_loss"),
    kernel = c("gaussian", "matern52", "matern32"))
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

  # Convert gamma to kernel hyperparameters (theta = 10^gamma)
  theta <- 10^gamma

  # Compute kernel matrix using specified kernel and theta
  K <- kernel_matrix(Xnorm, theta = theta, kernel = kernel)
  diag(K) <- 0  # Leave-One-Out Cross-Validation (LOOCV)

  if (model == "BKP") {
    ## -------- Binary case (Beta-Binomial) --------
    # get the prior parameters: alpha0(x) and beta0(x)
    prior_par <- get_prior(prior = prior, model = model, r0 = r0, p0 = p0, y = y, m = m, K = K)
    alpha0 <- prior_par$alpha0
    beta0 <- prior_par$beta0

    # Compute posterior alpha and beta
    alpha_n <- alpha0 + as.vector(K %*% y)
    beta_n <- beta0 + as.vector(K %*% (m - y))

    # Posterior mean prediction of success probability
    pi_hat <- alpha_n / (alpha_n + beta_n)
    pi_hat <- pmin(pmax(pi_hat, 1e-10), 1 - 1e-10)   # avoid log(0)

    # Empirical success rate
    pi_tilde <- y / m

    # ---------------- Loss computation ----------------
    if (loss == "brier") {
      # Standard Brier score (mean squared error)
      brier <- mean((pi_hat - pi_tilde)^2)
      return(brier)
    } else {
      # Standard log-loss (cross-entropy)
      # log_loss <- -mean(y * log(pi_hat) + (m - y) * log(1 - pi_hat))
      log_loss <- -mean(pi_tilde * log(pi_hat) + (1 - pi_tilde) * log(1 - pi_hat))
      return(log_loss)
    }
  } else {
    ## -------- Multiclass case (Dirichlet-Multinomial) --------
    # get the prior parameters: alpha0(x)
    alpha0 <- get_prior(prior = prior, model = model, r0 = r0, p0 = p0, Y = Y, K = K)

    # Compute posterior alpha
    alpha_n <- as.matrix(alpha0) + as.matrix(K %*% Y)

    # Posterior mean prediction of success probability
    pi_hat <- alpha_n / rowSums(alpha_n)
    pi_hat <- pmin(pmax(pi_hat, 1e-6), 1 - 1e-6)   # avoid log(0)

    # Empirical class probabilities
    pi_tilde <- Y / rowSums(Y)

    if (loss == "brier") {
      # Brier score (Mean Squared Error)
      brier <- mean((pi_hat - pi_tilde)^2)
      return(brier)
    } else {
      # log-loss (cross-entropy)
      log_loss <- -mean(pi_tilde * log(pi_hat))
      return(log_loss)
    }
  }
}
