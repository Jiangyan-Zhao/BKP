#' @title Loss Function for Fitting the DKP Model
#'
#' @description Computes the loss used to fit the DKP model. Supports the Brier
#'   score (mean squared error) and negative log-loss (cross-entropy), under
#'   different prior specifications.
#'
#' @inheritParams fit.DKP
#' @inheritParams loss_fun
#'
#' @return A single numeric value representing the total loss (to be minimized).
#'
#' @seealso \code{\link{loss_fun}}, \code{\link{fit.DKP}}, \code{\link{get_prior_dkp}},
#'   \code{\link{kernel_matrix}}
#'
#' @examples
#' set.seed(123)
#' n = 10
#' Xnorm = matrix(runif(2 * n), ncol = 2)
#' m = rep(10, n)
#' y = rbinom(n, size = m, prob = runif(n))
#' Y = cbind(y, m-y)
#' loss_fun_dkp(gamma = rep(0, 2), Xnorm = Xnorm, Y = Y)
#' @export


loss_fun_dkp <- function(
    gamma, Xnorm, Y,
    prior = c("noninformative", "fixed", "adaptive"), r0 = 2, p0 = NULL,
    loss = c("brier", "log_loss"),
    kernel = c("gaussian", "matern52", "matern32"))
{
  # Match and validate the loss and kernel arguments
  prior <- match.arg(prior)
  loss <- match.arg(loss)
  kernel <- match.arg(kernel)

  n <- nrow(Y)

  # Convert gamma to kernel hyperparameters (theta = 10^gamma)
  theta <- 10^gamma

  # Compute kernel matrix using specified kernel and theta
  K <- kernel_matrix(Xnorm, theta = theta, kernel = kernel)
  diag(K) <- 0  # Leave-One-Out Cross-Validation (LOOCV)

  # get the prior parameters: alpha0(x) and beta0(x)
  alpha0 <- get_prior_dkp(prior = prior, r0 = r0, p0 = p0, Y = Y, K = K)

  # Compute posterior alpha
  alpha_n <- as.matrix(alpha0) + as.matrix(K %*% Y)

  # Numerical stabilization: avoid log(0) or NaNs
  alpha_n <- pmax(alpha_n, 1e-10)

  # Posterior mean prediction of success probability
  pi_hat <- alpha_n / rowSums(alpha_n)
  pi_hat <- pmin(pmax(pi_hat, 1e-10), 1 - 1e-10)   # avoid log(0)

  if (loss == "brier") {
    # Brier score (Mean Squared Error)
    # Empirical success rate
    pi_tilde <- Y / rowSums(Y)
    # Brier score: mean squared error between predicted and observed
    brier <- mean((pi_hat - pi_tilde)^2)
    return(brier)
  } else if (loss == "log_loss"){
    # log-loss (cross-entropy)
    log_loss <- -mean( Y * log(pi_hat) )
    return(log_loss)
  } else {
    stop("Unsupported loss type. Use 'brier' or 'log_loss'.")
  }
}
