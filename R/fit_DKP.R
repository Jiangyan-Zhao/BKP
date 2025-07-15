#' @name fit.DKP
#'
#' @title Fit a Dirichlet Kernel Process (DKP) Model for Multinomial Observations
#'
#' @description
#' Fits a Dirichlet Kernel Process model for multinomiall data, estimating latent functions
#' representing Dirichlet distribution parameters using kernel smoothing.
#'
#' @param data Optional data frame containing covariates \code{X} and observations \code{Y}
#'   in the first \eqn{d} and \eqn{d+q} columns, respectively.
#' @param X Covariate matrix of size \eqn{n \times d}. Ignored if \code{data} is provided.
#' @param Y Matrix of observed successes of size \eqn{n \times q}.
#' @param Xbounds Optional \eqn{d \times 2} matrix defining lower and upper bounds for each input dimension
#'   (used for normalization). Defaults to \code{[0,1]^d}.
#' @param prior Type of prior used: "noninformative", "fixed", or "adaptive".
#' @param r0 Global precision parameter used in prior specification.
#' @param p0 Global prior mean used when prior = "fixed". Must be a vector of length \eqn{q}.
#' @param kernel Character string specifying the kernel function to use: one of \code{"gaussian"},
#'   \code{"matern52"}, or \code{"matern32"}.
#' @param loss Loss function to be minimized during hyperparameter tuning. Choose between \code{"brier"}
#'   (default) and \code{"NLML"} (negative log marginal likelihood).
#' @param num_multi_start Number of initial points for multi-start optimization (default = 10*d).
#'
#' @details
#' The function fits the DKP model and stores everything necessary for prediction.
#' The kernel hyper-parameters optimization uses multi-start gradient-based method.
#' Inputs are first normalized to the \eqn{[0,1]^d} space defined by \code{Xbounds}.
#' The kernel matrix is computed with optimized parameters, and posterior Dirichlet parameters
#' \eqn{\alpha_n} is computed by smoothing the observed successes \eqn{y} with the kernel matrix:
#' \deqn{
#' \alpha_n = \alpha_0 + K \cdot Y,
#' }
#' where \eqn{K} is the kernel matrix evaluated at normalized inputs.
#'
#' @return A list of class \code{"DKP"} containing the fitted model, with the following components:
#' \describe{
#'   \item{\code{bestTheta}}{A numeric vector of optimal kernel hyperparameters (lengthscale), obtained by minimizing the specified loss function.}
#'   \item{\code{kernel}}{Character string indicating the kernel type used ("gaussian", "matern52", or "matern32").}
#'   \item{\code{loss}}{Character string specifying the loss function used ("brier" or "log_loss").}
#'   \item{\code{minLoss}}{The minimum loss value achieved during hyperparameter optimization.}
#'   \item{\code{X}}{The original input matrix (unnormalized), with dimensions \code{n x d}.}
#'   \item{\code{Xnorm}}{The normalized input matrix scaled to \code{[0,1]^d}.}
#'   \item{\code{Xbounds}}{A \code{d x 2} matrix specifying the lower and upper bounds used for normalizing each dimension of \code{X}.}
#'   \item{\code{Y}}{A numeric matrix of observed successes (length \code{n x q}).}
#'   \item{\code{prior}}{Character string indicating the prior type used ("noninformative", "fixed", or "adaptive").}
#'   \item{\code{r0}}{Global prior precision parameter (only used for fixed or adaptive priors).}
#'   \item{\code{p0}}{Global prior mean parameter (only used for fixed priors).}
#'   \item{\code{alpha0}}{The prior alpha parameters for each input location. A vector or matrix with \code{n} rows.}
#'   \item{\code{alpha_n}}{Posterior alpha parameters computed as \code{alpha0 + K \%*\% Y}.}
#' }
#'
#' @author Jiangyan Zhao, Kunhai Qing, Jin Xu
#' @examples
#' ### 1D
#' set.seed(123)
#' n <- 30
#' Xbounds <- matrix(c(-2,2), nrow=1)
#' x <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' true_pi <- matrix(c(true_pi/2,true_pi/2,1-true_pi),nrow = n, byrow = F)
#' m <- sample(100, n, replace = TRUE)
#' Y <- matrix(0, nrow = n, ncol = 3)
#' for (i in 1:n) {
#'   Y[i, ] <- rmultinom(n=1, size=m[i], prob=true_pi[i, ])
#' }
#' DKPmodel <- fit.DKP(p=3,X=x,Y=Y,Xbounds = Xbounds,prior = "noninformative",kernel = "gaussian",loss = "brier")
#' print(DKPmodel)
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
#' true_pi <- matrix(c(true_pi/2,true_pi/2,1-true_pi),nrow = n, byrow = F)
#' m <- sample(100, n, replace = TRUE)
#' Y <- matrix(0, nrow = n, ncol = 3)
#' for (i in 1:n) {
#'   Y[i, ] <- rmultinom(n=1, size=m[i], prob=true_pi[i, ])
#' }
#' DKPmodel <- fit.DKP(p=3,X=x,Y=Y,Xbounds = Xbounds,prior = "noninformative",kernel = "gaussian",loss = "brier")
#' print(DKPmodel)
#'
#' @export
#' @importFrom tgp lhs
#' @importFrom optimx multistart

fit.DKP <- function(
    data = NULL, d = NULL, q = NULL, X = NULL, Y = NULL, Xbounds = NULL,
    prior = c("noninformative", "fixed", "adaptive"), r0 = NULL, p0 = NULL,
    kernel = c("gaussian", "matern52", "matern32"),
    loss = c("brier", "log_loss"),
    num_multi_start = NULL
){
  # Handle input data: prioritize 'data' data frame, otherwise use individual X and Y.
  if (!missing(data)) {
    if (ncol(data) < 3) {
      stop("The 'data' frame must contain at least three columns (covariates X and Y in order).")
    }
    if (is.null(q) & is.null(d)){
      stop("Either 'q' or 'd' must be provided.")
    }
    if (is.null(d)){
      d <- ncol(data) - q
    } else {
      if (is.null(q)) {
        q <- ncol(data) - d
      }
    }
    X <- as.matrix(data[, 1:d])
    Y <- as.matrix(data[, (d + 1):(d + q)])
  } else {
    if (is.null(X) || is.null(Y) ) {
      stop("Either 'data' must be provided, or both of 'X' and 'Y'.")
    }
    d <- ncol(X)
    q <- ncol(Y)
    X <- as.matrix(X)
    Y <- as.matrix(Y)
  }

  # Parse arguments
  prior <- match.arg(prior)
  kernel <- match.arg(kernel)
  loss <- match.arg(loss)

  # Fixed mean informative prior
  if (is.null(r0) || is.null(p0)){
    r0 <- q
    p0 <- as.vector(rep(1/q,q))
  }

  # Validity checks
  n <- nrow(X)
  if (nrow(Y) != n) {
    stop("The length of 'Y' must match the number of rows in X.")
  }
  # Normalize X to [0,1]^d
  if (is.null(Xbounds)) {
    Xbounds <- cbind(rep(0, d), rep(1, d))
  }
  if (!all(dim(Xbounds) == c(d, 2))) {
    stop("Xbounds must be a d x 2 matrix.")
  }
  Xnorm <- sweep(X, 2, Xbounds[,1], "-")
  Xnorm <- sweep(Xnorm, 2, Xbounds[,2] - Xbounds[,1], "/")

  # Generate initial gamma values using space-filling design followed by D-optimal selection.
  # Gamma corresponds to kernel scale via theta = 10^gamma.
  # Bounds are set based on input dimension d to ensure reasonable search range.
  gammaBounds <- matrix(c(rep((log10(d)+2)/2, d), rep((log10(d)-log10(500))/2, d)), ncol = 2)

  # Perform multi-start optimization to find the best kernel parameters.
  if(is.null(num_multi_start)){num_multi_start <- 10 * d}
  initialGamma <- tgp::lhs(num_multi_start, gammaBounds)
  res <- optimx::multistart(
    parmat = initialGamma,
    fn     = loss_fun_dkp,
    method = "L-BFGS-B",
    # lower  = gammaBounds[,1], upper  = gammaBounds[,2],
    lower  = rep(-10, d), upper  = rep(10, d),
    prior = prior, r0 = r0, p0 = p0,
    Xnorm = Xnorm, Y = Y,
    loss = loss, kernel = kernel,
    control= list(trace=0))

  # Extract the results from the optimization.
  bestIndex <- which.min(res$value) # Find the index of the minimum loss.
  bestGamma <- as.numeric(res[bestIndex, 1:d]) # Get the gamma parameters corresponding to min loss.
  bestTheta <- 10^(bestGamma) # Transform gamma back to the kernel parameters (theta).
  minLoss <- res$value[bestIndex] # Get the minimum loss value.

  # Compute kernel matrix with the optimized parameters
  K <- kernel_matrix(Xnorm, theta = bestTheta, kernel = kernel)

  # get the prior parameters: alpha0(x)
  prior_par <- get_prior_dkp(prior = prior, r0 = r0, p0 = p0, Y = Y, K = K)
  alpha0 <- prior_par$alpha0

  # Compute posterior Dirichlet parameters using kernel smoothing
  if (prior == "noninformative" || prior == "fixed"){
    alpha0 <- matrix(rep(alpha0,n),nrow = n , byrow = T)
  }
  alpha_n <- alpha0 + as.matrix(K %*% Y)

  # Construct the 'DKP' model object as a list.
  # This list contains all essential information about the fitted model.
  model <- list(
    bestTheta = bestTheta, kernel = kernel,
    loss = loss, minLoss = minLoss,
    X = X, Xnorm = Xnorm, Xbounds = Xbounds, Y = Y,
    prior = prior, r0 = r0, p0 = p0, alpha0 = alpha0,
    alpha_n = alpha_n
  )

  # Assign the "DKP" class to the object, enabling S3 generic methods.
  class(model) <- "DKP"
  return(model) # Return the fitted DKP model object.
}
