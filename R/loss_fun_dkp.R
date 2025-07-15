#' @name lossFun
#'
#' @title Loss Function for DKP Model Fitting
#'
#' @description
#' Computes the loss used for fitting the Beta Kernel Process (DKP) model.
#' Supports Brier score (mean squared error) or log-loss (cross-entropy),
#' under different prior specifications.
#'
#' @param gamma A numeric vector of kernel hyperparameters (on the log scale).
#' @param Xnorm A normalized input matrix (each column scaled \code{[0,1]}).
#' @param Y A matrix of observed successes of size \eqn{n \times q}.
#' @param prior Type of prior: "noninformative", "fixed", or "adaptive".
#' @param r0 Global prior strength (used in "fixed"/"adaptive").
#' @param p0 Prior mean for fixed prior.
#' @param loss Type of loss: "brier" or "log_loss".
#' @param kernel Type of kernel: "gaussian", "matern52", or "matern32".
#'
#' @return A scalar numeric value representing the loss.
#'
#' @author Jiangyan Zhao, Kunhai Qing, Jin Xu
#' @examples
#' loss_fun_dkp(gamma = rep(1, 2),
#'         Xnorm = matrix(runif(20), ncol=2),
#'         Y = matrix(runif(20), ncol=2),
#'         prior = "noninformative",
#'         loss = "brier",
#'         kernel = "gaussian")
#' @export


loss_fun_dkp <- function(gamma, Xnorm, Y,
                    prior = c("noninformative", "fixed", "adaptive"), r0, p0,
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
  prior_par <- get_prior_dkp(
    prior = prior, r0 = r0, p0 = p0,
    Y = Y, K = K
  )
  alpha0 <- prior_par$alpha0
  if (prior == "noninformative" || prior == "fixed"){
    alpha0 <- matrix(rep(alpha0,n),nrow = n , byrow = T)
  }
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
