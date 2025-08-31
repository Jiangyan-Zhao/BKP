#' @name predict
#'
#' @title Predict Method for BKP or DKP Models
#'
#' @description Generates posterior predictive summaries from a fitted BKP or
#'   DKP model at new input locations.
#'
#' @param object An object of class \code{"BKP"} or \code{"DKP"}, typically
#'   returned by \code{\link{fit.BKP}} or \code{\link{fit.DKP}}.
#' @param Xnew A numeric matrix (or vector) of new input locations where
#'   predictions are desired.
#' @param CI_level Credible level for prediction intervals (default is
#'   \code{0.95}, corresponding to 95% CI).
#' @param threshold Classification threshold for binary prediction based on
#'   posterior mean (used only for BKP; default is \code{0.5}).
#' @param ... Additional arguments passed to generic \code{predict} methods
#'   (currently not used; included for S3 method consistency).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{X}}{The training data points.}
#'   \item{\code{Xnew}}{The new data points.}
#'   \item{\code{mean}}{BKP: Posterior mean of the success probability at each location.
#'                      DKP: A matrix of posterior mean class probabilities (rows = inputs, columns = classes).}
#'   \item{\code{variance}}{BKP: Posterior variance of the success probability.
#'                          DKP: A matrix of posterior variances for each class.}
#'   \item{\code{lower}}{BKP: Lower bound of the prediction interval (e.g., 2.5th percentile for 95% CI).
#'                       DKP: A matrix of lower bounds for each class (e.g., 2.5th percentile).}
#'   \item{\code{upper}}{BKP: Upper bound of the prediction interval (e.g., 97.5th percentile for 95% CI).
#'                       DKP: A matrix of upper bounds for each class (e.g., 97.5th percentile).}
#'   \item{\code{class}}{BKP: Predicted binary label (0 or 1), based on posterior mean and threshold, if \code{m = 1}.
#'                       DKP: Predicted class label (i.e., the class with the highest posterior mean probability).}
#' }
#'
#' @seealso \code{\link{fit.BKP}} for fitting Beta Kernel Process models.
#'   \code{\link{fit.DKP}} for fitting Dirichlet Kernel Process models.
#'   \code{\link{plot.BKP}} for visualizing fitted BKP models.
#'   \code{\link{plot.DKP}} for visualizing fitted DKP models.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv.
#'   https://doi.org/10.48550/arXiv.2508.10447.
#'
#' @keywords BKP
#'
#' @examples
#' # ============================================================== #
#' # ========================= BKP Examples ======================= #
#' # ============================================================== #
#'
#' #-------------------------- 1D Example ---------------------------
#' set.seed(123)
#'
#' # Define true success probability function
#' true_pi_fun <- function(x) {
#'   (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' }
#'
#' n <- 30
#' Xbounds <- matrix(c(-2,2), nrow=1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit BKP model
#' model1 <- fit.BKP(X, y, m, Xbounds=Xbounds)
#'
#' # Prediction on training data
#' predict(model1)
#'
#' # Prediction on new data
#' Xnew = matrix(seq(-2, 2, length = 10), ncol=1) #new data points
#' predict(model1, Xnew = Xnew)
#'
#'
#' #-------------------------- 2D Example ---------------------------
#' set.seed(123)
#'
#' # Define 2D latent function and probability transformation
#' true_pi_fun <- function(X) {
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
#'   return(pnorm(f))  # Transform to probability
#' }
#'
#' n <- 100
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit BKP model
#' model2 <- fit.BKP(X, y, m, Xbounds=Xbounds)
#'
#' # Prediction on training data
#' predict(model2)
#'
#' # Prediction on new data
#' x1 <- seq(Xbounds[1,1], Xbounds[1,2], length.out = 10)
#' x2 <- seq(Xbounds[2,1], Xbounds[2,2], length.out = 10)
#' Xnew <- expand.grid(x1 = x1, x2 = x2)
#' predict(model2, Xnew = Xnew)
#'
#' @export
#' @method predict BKP

predict.BKP <- function(object, Xnew = NULL, CI_level = 0.95, threshold = 0.5, ...)
{
  # Extract components
  X       <- object$X
  Xnorm   <- object$Xnorm
  y       <- object$y
  m       <- object$m
  theta   <- object$theta_opt
  kernel  <- object$kernel
  prior   <- object$prior
  r0      <- object$r0
  p0      <- object$p0
  Xbounds <- object$Xbounds
  d       <- ncol(Xnorm)

  if(!is.null(Xnew)){
    # Ensure Xnew is a matrix and matches input dimension
    if (is.null(nrow(Xnew))) {
      Xnew <- matrix(Xnew, nrow = 1)
    }
    Xnew <- as.matrix(Xnew)
    if (ncol(Xnew) != d) {
      stop("The number of columns in 'Xnew' must match the original input dimension.")
    }

    # Normalize Xnew to [0,1]^d
    Xnew_norm <- sweep(Xnew, 2, Xbounds[, 1], "-")
    Xnew_norm <- sweep(Xnew_norm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

    # Compute kernel matrix
    K <- kernel_matrix(Xnew_norm, Xnorm, theta = theta, kernel = kernel) # m*n matrix
  }else{
    # Compute kernel matrix
    K <- kernel_matrix(Xnorm, theta = theta, kernel = kernel) # n*n matrix
  }

  # get the prior parameters: alpha0(x) and beta0(x)
  prior_par <- get_prior(prior = prior, model = "BKP",
                         r0 = r0, p0 = p0, y = y, m = m, K = K)
  alpha0 <- prior_par$alpha0
  beta0 <- prior_par$beta0

  # Posterior parameters
  alpha_n <- pmax(alpha0 + as.vector(K %*% y), 1e-6)
  beta_n  <- pmax(beta0 + as.vector(K %*% (m - y)), 1e-6)

  # Predictive mean and variance
  pi_mean <- alpha_n / (alpha_n + beta_n)
  pi_var  <- pi_mean * (1 - pi_mean) / (alpha_n + beta_n + 1)

  # Credible intervals
  pi_lower <- qbeta((1 - CI_level) / 2, alpha_n, beta_n)
  pi_upper <- qbeta((1 + CI_level) / 2, alpha_n, beta_n)


  # Output list
  prediction <- list(
    X        = X,
    Xnew     = Xnew,
    mean     = pi_mean,
    variance = pi_var,
    lower    = pi_lower,
    upper    = pi_upper,
    CI_level  = CI_level)

  # Posterior classification label (only for classification data)
  if (all(m == 1)) {
    prediction$class <- ifelse(pi_mean > threshold, 1, 0)
    prediction$threshold <- threshold
  }

  class(prediction) <- "predict.BKP"
  return(prediction)
}
