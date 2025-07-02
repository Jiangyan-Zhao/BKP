#' @name fit.BKP
#'
#' @title Fit a Beta Kernel Process (BKP) Model for Binomial Observations
#'
#' @description
#' Fits a Beta Kernel Process model for binomial data, estimating latent functions
#' representing Beta distribution parameters using kernel smoothing.
#'
#' @param data Optional data frame containing covariates \code{X}, observations \code{y}, and counts \code{m}
#'   in the first \eqn{d}, \eqn{d+1}, and \eqn{d+2} columns, respectively.
#' @param X Covariate matrix of size \eqn{n \times d}. Ignored if \code{data} is provided.
#' @param y Vector of observed successes. Must be numeric with length equal to number of rows in \code{X}.
#' @param m Vector of binomial counts (trials). Must be numeric, positive, and same length as \code{y}.
#' @param alpha0 Prior shape parameter (scalar or vector of length \eqn{n}) for the beta distribution.
#' @param beta0 Prior shape parameter (scalar or vector of length \eqn{n}) for the beta distribution.
#' @param Xbounds Optional \eqn{d \times 2} matrix defining lower and upper bounds for each input dimension
#'   (used for normalization). Defaults to \code{[0,1]^d}.
#' @param kernel Character string specifying the kernel function to use: one of \code{"gaussian"},
#'   \code{"matern52"}, or \code{"matern32"}.
#' @param loss Loss function to be minimized during hyperparameter tuning. Choose between \code{"brier"}
#'   (default) and \code{"NLML"} (negative log marginal likelihood).
#'
#' @details
#' The function fits the BKP model by optimizing kernel hyperparameters using multi-start
#' local optimization. Inputs are first normalized to the \eqn{[0,1]^d} space defined by \code{Xbounds}.
#' The kernel matrix is computed with optimized parameters, and posterior Beta parameters
#' \eqn{\alpha_n} and \eqn{\beta_n} are computed by smoothing the observed successes \eqn{y}
#' and failures \eqn{m - y} with the kernel matrix:
#' \deqn{
#' \alpha_n = \alpha_0 + K \cdot y, \quad
#' \beta_n = \beta_0 + K \cdot (m - y)
#' }
#' where \eqn{K} is the kernel matrix evaluated at normalized inputs.
#'
#' @return An object of class \code{"BKP"} containing:
#'   \item{bestTheta}{Optimal kernel parameters (length \eqn{d}).}
#'   \item{kernel}{Kernel type used.}
#'   \item{loss}{Loss function used.}
#'   \item{minLoss}{Minimum achieved loss.}
#'   \item{X}{Original (unnormalized) design matrix.}
#'   \item{Xnorm}{Normalized design matrix (scaled to \code{[0,1]^d}).}
#'   \item{Xbounds}{\eqn{d \times 2} matrix defining lower and upper bounds for each input dimension.}
#'   \item{y}{Observed successes.}
#'   \item{m}{Binomial counts.}
#'   \item{alpha0, beta0}{Prior parameters used.}
#'   \item{alpha.n, beta.n}{Posterior parameters (vector).}
#'
#' @author Jiangyan Zhao, Kunhai Qing, Jin Xu
#'
#' @examples
#' ### 1D
#' set.seed(123)
#' n <- 30
#' Xbounds <- matrix(c(-2,2), nrow=1)
#' x <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#' df <- data.frame(x = x, y = y, m = m)
#' xx = matrix(seq(-2, 2, length = 100), ncol=1) # new data points
#' model <- fit.BKP(df, Xbounds=Xbounds)
#' print(model)
#'
#' ### 2D
#' set.seed(123)
#' n <- 100
#' f <- function(X) {
#'   if(is.null(nrow(X))) X <- matrix(X, nrow=1)
#'   m <- 8.6928
#'   s <- 2.4269
#'   x1 <- 4*X[,1]- 2
#'   x2 <- 4*X[,2]- 2
#'   a <- 1 + (x1 + x2 + 1)^2 *
#'     (19- 14*x1 + 3*x1^2- 14*x2 + 6*x1*x2 + 3*x2^2)
#'   b <- 30 + (2*x1- 3*x2)^2 *
#'     (18- 32*x1 + 12*x1^2 + 48*x2- 36*x1*x2 + 27*x2^2)
#'   f <- log(a*b)
#'   f <- (f- m)/s
#'   return(f) }
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' x <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- pnorm(f(x))
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#' df <- data.frame(x = x, y = y, m = m)
#' xx1 <- seq(Xbounds[1,1], Xbounds[1,2], length.out = 100)
#' xx2 <- seq(Xbounds[2,1], Xbounds[2,2], length.out = 100)
#' xx <- expand.grid(xx1 = xx1, xx2 = xx2)
#' #plot the true probability
#' true_pi <- matrix(pnorm(f(xx)), nrow = length(xx1), ncol = length(xx2))
#' image(xx1, xx2, true_pi, xlab ="X1", ylab ="X2",
#'                 main = "True Probability",
#'                 col = hcl.colors(100, "viridis"))
#' contour(xx1, xx2, true_pi, add = TRUE, col = "black")
#' model <- fit.BKP(df)
#' print(model)
#'
#' @export
#' @importFrom tgp lhs dopt.gp
#' @importFrom optimx multistart

fit.BKP <- function(
    data, X = NULL, y = NULL, m = NULL,
    alpha0 = 1, beta0 = 1,
    Xbounds = NULL,
    kernel = c("gaussian", "matern52", "matern32"),
    loss = c("brier", "log_Loss", "NLML")
){

  # Handle input data: prioritize 'data' data frame, otherwise use individual X, y, and m.
  if (!missing(data)) {
    if (ncol(data) < 3) {
      stop("The 'data' frame must contain at least three columns (covariates X, y, and m in order).")
    }
    d <- ncol(data) - 2
    X <- as.matrix(data[, 1:d])
    y <- as.matrix(data[, d + 1])
    m <- as.matrix(data[, d + 2])
  } else {
    if (is.null(X) || is.null(y) || is.null(m)) {
      stop("Either the 'data' argument must be supplied, or 'X', 'y', and 'm' must all be provided separately.")
    }

    X <- as.matrix(X)
    d <- ncol(X)
    y <- matrix(y, ncol = 1)
    m <- matrix(m, ncol = 1)
  }

  # Type checks
  if (!all(sapply(X, is.numeric))) stop("All columns of X (covariates) must be numeric.")
  if (!all(sapply(y, is.numeric))) stop("'y' must be numeric.")
  if (!all(sapply(m, is.numeric))) stop("'m' must be numeric.")

  # Dimension consistency checks
  n <- nrow(X)
  if (nrow(y) != n || nrow(m) != n) {
    stop("The lengths of 'y' and 'm' must match the number of rows in X.")
  }

  # Value validity checks
  if (any(y < 0) || any(m <= 0) || any(y > m)) {
    stop("Each value of 'y' must be between 0 and its corresponding 'm', and 'm' must be strictly positive.")
  }

  # Default Xbounds to [0,1]^d if not provided
  if (is.null(Xbounds)) {
    Xbounds <- matrix(c(rep(0, d), rep(1, d)), ncol = 2)
  }

  if (nrow(Xbounds) != d || ncol(Xbounds) != 2) {
    stop("Xbounds must be a matrix with dimension d x 2 (rows = number of features)")
  }

  # # Ensure alpha0 and beta0 are either scalars or vectors of length n; expand scalars if needed
  # if (length(alpha0) == 1) {
  #   alpha0 <- rep(alpha0, n)
  # } else if (length(alpha0) != n) {
  #   stop("alpha0 must be either a scalar or a vector of length equal to the number of observations.")
  # }
  #
  # if (length(beta0) == 1) {
  #   beta0 <- rep(beta0, n)
  # } else if (length(beta0) != n) {
  #   stop("beta0 must be either a scalar or a vector of length equal to the number of observations.")
  # }


  # Ensure alpha0 and beta0 are either scalars or vectors of length n
  if (length(alpha0) != 1 && length(alpha0) != n) {
    stop("alpha0 must be either a scalar or a vector of length equal to the number of observations.")
  }

  if (length(beta0) != 1 && length(beta0) != n) {
    stop("beta0 must be either a scalar or a vector of length equal to the number of observations.")
  }

  # Match the kernel argument explicitly
  kernel <- match.arg(kernel)

  # Match the loss argument explicitly
  loss <- match.arg(loss)

  # Normalize X to [0, 1]^d
  Xnorm <- sweep(X, 2, Xbounds[,1], "-")
  Xnorm <- sweep(Xnorm, 2, Xbounds[,2] - Xbounds[,1], "/")

  # Generate initial gamma values using space-filling design followed by D-optimal selection.
  # Gamma corresponds to kernel scale via theta = 10^gamma.
  # Bounds are set based on input dimension d to ensure reasonable search range.
  gammaBounds <- matrix(
    c(rep((log10(d)+2)/2, d), rep((log10(d)-log10(500))/2, d)),
    ncol = 2
  )
  # Latin Hypercube Sampling in gamma space
  gammaCandidates <- lhs(100 * d, gammaBounds)
  # Select 10*d initial values using sequential D-optimal design
  initialGamma <- dopt.gp(10 * d, Xcand = gammaCandidates)$XX

  # Perform multi-start optimization to find the best kernel parameters.
  res <- multistart(
    parmat = initialGamma,
    fn     = lossFun,
    method = "L-BFGS-B",
    # lower  = gammaBounds[,1], upper  = gammaBounds[,2],
    lower  = rep(-10, d), upper  = rep(10, d),
    alpha0 = alpha0, beta0 = beta0,
    Xnorm = Xnorm, y = y, m=m,
    loss = loss, kernel = kernel,
    control= list(trace=0))

  # Extract the results from the optimization.
  bestIndex <- which.min(res$value) # Find the index of the minimum loss.
  bestGamma <- as.numeric(res[bestIndex, 1:d]) # Get the gamma parameters corresponding to min loss.
  bestTheta <- 10^(bestGamma) # Transform gamma back to the kernel parameters (theta).
  minLoss <- res$value[bestIndex] # Get the minimum loss value.

  # Compute kernel matrix with the optimized parameters
  K <- kernel_matrix(Xnorm, theta = bestTheta, kernel = kernel)
  # Compute posterior Beta parameters using kernel smoothing
  alpha.n <- alpha0 + as.vector(t(K) %*% y)
  beta.n <- beta0 + as.vector(t(K) %*% (m - y))

  # Construct the 'BKP' model object as a list.
  # This list contains all essential information about the fitted model.
  BKPmodel <- list(bestTheta = bestTheta, kernel = kernel,
    loss = loss, minLoss = minLoss,
    X = X, Xnorm = Xnorm, Xbounds = Xbounds, y = y, m = m,
    alpha0 = alpha0, beta0 = beta0, alpha.n = alpha.n, beta.n = beta.n
  )

  # Assign the "BKP" class to the object, enabling S3 generic methods.
  class(BKPmodel) <- "BKP"
  return(BKPmodel) # Return the fitted BKP model object.
}
