#' @name lossFun
#'
#' @title Loss Function for BKP Model Fitting
#'
#' @description
#' Computes the loss used for fitting the Beta Kernel Process (BKP) model.
#' Supports Brier score or negative log marginal likelihood (NLML).
#'
#' @param gamma A numeric vector of kernel hyperparameters (on the log scale).
#' @param Xnorm A normalized input matrix (each column scaled to the range 0-1).
#' @param y A numeric vector of binary counts (successes).
#' @param m A numeric vector of total trials (denominators).
#' @param alpha0 Prior alpha parameter(s) for the Beta prior. Can be a scalar or a vector of length \code{length(y)}.
#' @param beta0 Prior beta parameter(s) for the Beta prior. Can be a scalar or a vector of length \code{length(y)}.
#' @param loss Character string indicating the loss function: either \code{"brier"} or \code{"NLML"}.
#' @param kernel Type of kernel function used in the BKP model. One of \code{"gaussian"}, \code{"matern52"}, or \code{"matern32"}.
#'
#' @return A scalar numeric value representing the loss.
#'
#' @author Jiangyan Zhao, Kunhai Qing, Jin Xu
#'
#' @examples
#' lossFun(gamma = rep(0, 2), Xnorm = matrix(runif(20), ncol=2),
#'         y = rbinom(10, 10, 0.5), m = rep(10, 10),
#'         alpha0 = 1, beta0 = 1, loss = "brier", kernel = "gaussian")
#'
#' @export

lossFun <- function(gamma, Xnorm, y, m, alpha0 = 1, beta0 = 1,
                    loss = c("brier", "NLML"),
                    kernel = c("gaussian", "matern52", "matern32"))
{
  # Match and validate the loss and kernel arguments
  loss <- match.arg(loss)
  kernel <- match.arg(kernel)

  n <- length(y)

  # Expand scalar alpha0 and beta0 to vectors of length n
  if (length(alpha0) == 1) alpha0 <- rep(alpha0, n)
  if (length(beta0)  == 1) beta0  <- rep(beta0, n)

  # Convert gamma to kernel hyperparameters (theta = 10^(-gamma))
  theta <- 10^(-gamma)

  # Compute kernel matrix using specified kernel and theta
  K <- kernel_matrix(Xnorm, theta = theta, kernel = kernel)

  # Compute posterior alpha and beta
  alpha.n <- alpha0 + as.vector(K %*% y)
  beta.n <- beta0 + as.vector(K %*% (m - y))

  # Numerical stabilization: avoid log(0) or NaNs
  alpha.n <- pmax(alpha.n, 1e-8)
  beta.n  <- pmax(beta.n, 1e-8)

  if (loss == "brier") {
    # Posterior mean prediction of success probability
    pi.tilde <- (alpha.n - y) / (alpha.n + beta.n - m)
    # Empirical success rate
    pi.hat <- y / m
    # Brier score: mean squared error between predicted and observed
    brier <- mean((pi.tilde - pi.hat)^2)
    return(brier)

  } else if (loss == "NLML") {
    # Negative log marginal likelihood (based on beta-binomial conjugacy)
    nlml <- sum(
      lbeta(alpha.n + y, beta.n + m - y) -
        lbeta(alpha.n, beta.n) +
        lchoose(m, y)
    )
    return(-nlml)  # Minimization target
  }
}
