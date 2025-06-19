#' @name predict
#'
#' @title Predict method for Beta Kernel Process (BKP) models
#'
#' @description
#' Generate predictions from a fitted BKP model at new input locations.
#'
#' @param object A fitted BKP model object returned by \code{\link{fit.BKP}}.
#' @param Xnew A matrix (or vector) of new input points at which to predict.
#' @param CI.size Confidence interval level (default = 0.05 for 95% CI).
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{X}{Original new input locations}
#'   \item{mean}{Posterior mean of success probability}
#'   \item{variance}{Posterior variance}
#'   \item{lower}{Lower bound of CI for success probability}
#'   \item{upper}{Upper bound of CI for success probability}
#' }
#'
#' @author Jiangyan Zhao, Kunhai Qing, Jin Xu
#'
#' @keywords BKP
#'
#' @export
#' @method predict BKP


predict.BKP <- function(object, Xnew, CI.size = 0.05)
{
  if (!inherits(object, "BKP")) {
    stop("The input is not of class 'BKP'. Please provide a model fitted with 'fit.BKP()'.")
  }

  BKPmodel <- object

  # Extract components
  Xnorm   <- BKPmodel$Xnorm
  y       <- BKPmodel$y
  m       <- BKPmodel$m
  theta   <- BKPmodel$bestTheta
  kernel  <- BKPmodel$kernel
  alpha0  <- BKPmodel$alpha0
  beta0   <- BKPmodel$beta0
  Xbounds <- BKPmodel$Xbounds
  d       <- ncol(Xnorm)

  # Ensure Xnew is a matrix and matches input dimension
  if (is.null(nrow(Xnew))) {
    Xnew <- matrix(Xnew, nrow = 1)
  }
  Xnew <- as.matrix(Xnew)
  if (ncol(Xnew) != d) {
    stop("The number of columns in 'Xnew' must match the original input dimension.")
  }

  # Normalize Xnew to [0,1]^d
  XnewNorm <- sweep(Xnew, 2, Xbounds[, 1], "-")
  XnewNorm <- sweep(XnewNorm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

  # Compute kernel matrix
  K <- kernel_matrix(Xnorm, XnewNorm, theta = theta, kernel = kernel)

  # Posterior parameters
  alpha.n <- alpha0 + as.vector(t(K) %*% y)
  beta.n  <- beta0 + as.vector(t(K) %*% (m - y))

  # Predictive mean and variance
  pi.mean <- alpha.n / (alpha.n + beta.n)
  pi.var  <- pi.mean * (1 - pi.mean) / (alpha.n + beta.n + 1)

  # Confidence intervals
  pi.lower <- qbeta(CI.size / 2, alpha.n, beta.n)
  pi.upper <- qbeta(1 - CI.size / 2, alpha.n, beta.n)

  prediction <- data.frame(X = I(Xnew), mean = pi.mean, variance = pi.var,
                           lower = pi.lower, upper = pi.upper)

  return(prediction)
}
