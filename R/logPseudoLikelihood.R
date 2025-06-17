#' @title Calculate log pseudo likelihood for BKP Model Optimization
#'
#' @description
#' This internal function computes minus log pseudo likelihood, which serves as the loss function
#' during the optimization process of the Beta Kernel Process (BKP) model.
#'
#' @param gamma A numeric vector representing the log-transformed kernel parameters.
#'   These are converted back to the actual kernel parameters (`theta`) using `10^(-gamma)`
#'   before being passed to the kernel function.
#' @param alpha0 A scalar or vector specifying the prior parameter for the beta distribution (shape1).
#' @param beta0 A scalar or vector specifying the prior parameter for the beta distribution (shape2).
#' @param kfun A kernel function previously defined within `fit.bkp()`. This function
#'   calculates the kernel similarity between two data points given kernel parameters.
#' @param X A matrix or data frame of covariate data, typically the normalized `X` from `fit.bkp()`.
#' @param y A numeric vector representing the number of successes for each observation.
#' @param m A numeric vector representing the total number of trials for each observation.
#' @param n An integer representing the sample size (number of observations).
#'
#' @return A single numeric value representing the calculated log pseudo likelihood.
#'   This score is minimized during the BKP model fitting process.
#'
#' @seealso
#' \code{\link{fit.bkp}} which uses this function for optimization.
#'
#' @keywords internal
#' This function is not exported for direct user use.
#'

logPseudoLikelihood <- function(gamma, alpha0, beta0, kfun, X, y, m, n){
  # Define a wrapper function for the kernel function.
  # This wrapper allows `kfun` to be used with `outer()` by transforming 'gamma'
  # back to 'theta' (kernel parameters) as `10^(-gamma)`.
  kernelWrapper <- function(i, j) {
    kfun(X[i, ], X[j, ], 10^(-gamma))
  }

  # Calculate the full N x N kernel matrix (K).
  # `outer()` applies `kernelWrapper` for all pairs of indices (i, j) from 1 to n.
  # `Vectorize()` ensures that `kernelWrapper` works correctly with `outer()`
  # even if it's not inherently vectorized for pair-wise operations.
  K <- base::outer(
    1:n, 1:n,
    base::Vectorize(kernelWrapper)
  )
  # Calculate the updated alpha (Alpha) and alpha + beta (AlphaPlusBeta) parameters
  # for the Beta-Binomial predictive distribution.
  # These updates incorporate the prior beliefs (alpha0, beta0) with observed data (y, m-y)
  # weighted by the kernel similarities (K).
  Alpha <- alpha0 + K %*% y
  Beta <- beta0 + K %*% (m-y)

  # Compute the Brier score.
  # This is the sum of squared differences between the predicted probabilities
  # (Alpha / AlphaPlusBeta, which is the mean of the posterior Beta distribution)
  # and the observed proportions (y / m).
  -base::sum( lbeta(Alpha+y,Beta+m-y) - lbeta(Alpha,Beta) )
}
