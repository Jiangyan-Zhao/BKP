#' @title Loss Function for Fitting the BKP Model
#'
#' @description Computes the loss used to fit the BKP model. Supports the Brier
#'   score (mean squared error) and negative log-loss (cross-entropy), under
#'   different prior specifications.
#'
#' @inheritParams fit.BKP
#' @param gamma A numeric vector of log-transformed kernel hyperparameters.
#' @param Xnorm A numeric matrix of normalized inputs (each column scaled to
#'   \code{[0,1]}).
#'
#' @return A single numeric value representing the total loss (to be minimized).
#'
#' @seealso \code{\link{loss_fun_dkp}}, \code{\link{fit.BKP}}, \code{\link{get_prior}},
#'   \code{\link{kernel_matrix}}
#'
#' @examples
#' set.seed(123)
#' n = 10
#' Xnorm = matrix(runif(2 * n), ncol = 2)
#' m = rep(10, n)
#' y = rbinom(n, size = m, prob = runif(n))
#' loss_fun(gamma = rep(0, 2), Xnorm = Xnorm, y = y, m = m)
#'
#' @export

loss_fun <- function(
    gamma, Xnorm, y, m,
    prior = c("noninformative", "fixed", "adaptive"), r0 = 2, p0 = 0.5,
    loss = c("brier", "log_loss"),
    kernel = c("gaussian", "matern52", "matern32"))
{
  # Match and validate the loss and kernel arguments
  prior <- match.arg(prior)
  loss <- match.arg(loss)
  kernel <- match.arg(kernel)

  n <- length(y)

  # Convert gamma to kernel hyperparameters (theta = 10^gamma)
  theta <- 10^gamma

  # Compute kernel matrix using specified kernel and theta
  K <- kernel_matrix(Xnorm, theta = theta, kernel = kernel)
  diag(K) <- 0  # Leave-One-Out Cross-Validation (LOOCV)

  # get the prior parameters: alpha0(x) and beta0(x)
  prior_par <- get_prior(
    prior = prior, r0 = r0, p0 = p0,
    y = y, m = m, K = K
  )
  alpha0 <- prior_par$alpha0
  beta0 <- prior_par$beta0

  # Compute posterior alpha and beta
  alpha_n <- alpha0 + as.vector(K %*% y)
  beta_n <- beta0 + as.vector(K %*% (m - y))

  # Numerical stabilization: avoid log(0) or NaNs
  alpha_n <- pmax(alpha_n, 1e-10)
  beta_n  <- pmax(beta_n, 1e-10)

  # Posterior mean prediction of success probability
  pi_hat <- alpha_n / (alpha_n + beta_n)
  pi_hat <- pmin(pmax(pi_hat, 1e-10), 1 - 1e-10)   # avoid log(0)

  if (loss == "brier") {
    # Brier score (Mean Squared Error)
    # Empirical success rate
    pi_tilde <- y / m
    # Brier score: mean squared error between predicted and observed
    brier <- mean((pi_hat - pi_tilde)^2)
    return(brier)
  } else if (loss == "log_loss"){
    # log-loss (cross-entropy)
    log_loss <- -mean(y * log(pi_hat) + (m - y) * log(1 - pi_hat))
    return(log_loss)
  } else {
    stop("Unsupported loss type. Use 'brier' or 'log_loss'.")
  }
}
