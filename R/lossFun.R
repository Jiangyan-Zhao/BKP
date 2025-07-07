#' @name lossFun
#'
#' @title Loss Function for BKP Model Fitting
#'
#' @description
#' Computes the loss used for fitting the Beta Kernel Process (BKP) model.
#' Supports Brier score (mean squared error) or log-loss (cross-entropy),
#' under different prior specifications.
#'
#' @param gamma A numeric vector of kernel hyperparameters (on the log scale).
#' @param Xnorm A normalized input matrix (each column scaled \code{[0,1]}).
#' @param y A numeric vector of binary counts (successes).
#' @param m A numeric vector of total trials (denominators).
#' @param prior Type of prior: "noninformative", "fixed", or "adaptive".
#' @param r0 Global prior strength (used in "fixed"/"adaptive").
#' @param p0 Prior mean for fixed prior (in (0,1)).
#' @param loss Type of loss: "brier" or "log_loss".
#' @param kernel Type of kernel: "gaussian", "matern52", or "matern32".
#'
#' @return A scalar numeric value representing the loss.
#'
#' @author Jiangyan Zhao, Kunhai Qing, Jin Xu
#'
#' @examples
#' lossFun(gamma = rep(0, 2),
#'         Xnorm = matrix(runif(20), ncol=2),
#'         y = rbinom(10, 10, 0.5),
#'         m = rep(10, 10))
#'
#' @export

lossFun <- function(gamma, Xnorm, y, m,
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
