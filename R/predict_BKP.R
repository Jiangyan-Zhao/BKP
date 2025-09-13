#' @name predict
#'
#' @title Posterior Prediction for BKP or DKP Models
#'
#' @description Generates posterior predictive summaries from a fitted Beta
#'   Kernel Process (BKP) or Dirichlet Kernel Process (DKP) model at new input
#'   locations. Supports prediction of posterior mean, variance, credible
#'   intervals, and classification labels (where applicable).
#'
#' @param object An object of class \code{"BKP"} or \code{"DKP"}, typically
#'   returned by \code{\link{fit_BKP}} or \code{\link{fit_DKP}}.
#' @param Xnew A numeric vector or matrix of new input locations at which to
#'   generate predictions. If \code{NULL}, predictions are returned for the
#'   training data.
#' @param CI_level Numeric between 0 and 1 specifying the credible level for
#'   posterior intervals (default \code{0.95} for 95% credible interval).
#' @param threshold Numeric between 0 and 1 specifying the classification
#'   threshold for binary predictions based on posterior mean (used only for
#'   BKP; default is \code{0.5}).
#' @param ... Additional arguments passed to generic \code{predict} methods
#'   (currently not used; included for S3 method consistency).
#'
#' @return A list containing posterior predictive summaries:
#' \describe{
#'   \item{\code{X}}{The original training input locations.}
#'   \item{\code{Xnew}}{The new input locations for prediction (same as \code{Xnew} if provided).}
#'   \item{\code{alpha_n}, \code{beta_n}}{Posterior shape parameters for each location:
#'     \itemize{
#'       \item BKP: Vectors of posterior shape parameters (\code{alpha_n}, \code{beta_n}) for each input location.
#'       \item DKP: Posterior shape parameter matrix \code{alpha_n} (rows = input locations, columns = classes).
#'     }}
#'   \item{\code{mean}}{Posterior mean prediction:
#'     \itemize{
#'       \item BKP: Posterior mean success probability at each location.
#'       \item DKP: Matrix of posterior mean class probabilities (rows = inputs, columns = classes).
#'     }}
#'   \item{\code{variance}}{Posterior predictive variance:
#'     \itemize{
#'       \item BKP: Variance of success probability.
#'       \item DKP: Matrix of predictive variances for each class.
#'     }}
#'   \item{\code{lower}}{Lower bound of the posterior credible interval:
#'     \itemize{
#'       \item BKP: Lower bound (e.g., 2.5th percentile for 95% CI).
#'       \item DKP: Matrix of lower bounds for each class.
#'     }}
#'   \item{\code{upper}}{Upper bound of the posterior credible interval:
#'     \itemize{
#'       \item BKP: Upper bound (e.g., 97.5th percentile for 95% CI).
#'       \item DKP: Matrix of upper bounds for each class.
#'     }}
#'   \item{\code{class}}{Predicted label:
#'     \itemize{
#'       \item BKP: Binary class (0 or 1) based on posterior mean and \code{threshold}, only if \code{m = 1}.
#'       \item DKP: Predicted class label with highest posterior mean probability.
#'     }}
#'   \item{\code{CI_level}}{The specified credible interval level.}
#' }
#'
#' @seealso \code{\link{fit_BKP}} and \code{\link{fit_DKP}} for model fitting;
#'   \code{\link{plot.BKP}} and \code{\link{plot.DKP}} for visualization of
#'   fitted models.
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
#' model1 <- fit_BKP(X, y, m, Xbounds=Xbounds)
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
#' model2 <- fit_BKP(X, y, m, Xbounds=Xbounds)
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
  #---- Argument checks ----
  X       <- object$X
  d       <- ncol(X)
  if (!is.null(Xnew)) {
    if (is.null(nrow(Xnew))) {
      Xnew <- matrix(Xnew, nrow = 1)
    }
    Xnew <- as.matrix(Xnew)
    if (!is.numeric(Xnew)) {
      stop("'Xnew' must be numeric.")
    }
    if (ncol(Xnew) != d) {
      stop("The number of columns in 'Xnew' must match the original input dimension.")
    }
  }
  if (!is.numeric(CI_level) || length(CI_level) != 1 || CI_level <= 0 || CI_level >= 1) {
    stop("'CI_level' must be a single numeric value strictly between 0 and 1.")
  }
  if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0 || threshold >= 1) {
    stop("'threshold' must be a single numeric value strictly between 0 and 1.")
  }


  if(!is.null(Xnew)){
    # Extract components
    Xnorm   <- object$Xnorm
    y       <- object$y
    m       <- object$m
    theta   <- object$theta_opt
    kernel  <- object$kernel
    prior   <- object$prior
    r0      <- object$r0
    p0      <- object$p0
    Xbounds <- object$Xbounds

    # Normalize Xnew to [0,1]^d
    Xnew_norm <- sweep(Xnew, 2, Xbounds[, 1], "-")
    Xnew_norm <- sweep(Xnew_norm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

    # Compute kernel matrix
    K <- kernel_matrix(Xnew_norm, Xnorm, theta = theta, kernel = kernel) # m*n matrix

    # get the prior parameters: alpha0(x) and beta0(x)
    prior_par <- get_prior(prior = prior, model = "BKP",
                           r0 = r0, p0 = p0, y = y, m = m, K = K)
    alpha0 <- prior_par$alpha0
    beta0 <- prior_par$beta0

    # Posterior parameters
    alpha_n <- alpha0 + as.vector(K %*% y)
    beta_n  <- beta0 + as.vector(K %*% (m - y))
  }else{
    # Use training data
    alpha_n <- object$alpha_n
    beta_n  <- object$beta_n
    m <- object$m
  }

  # Predictive mean and variance
  pred_mean <- alpha_n / (alpha_n + beta_n)
  pred_var  <- pred_mean * (1 - pred_mean) / (alpha_n + beta_n + 1)

  # Credible intervals
  pred_lower <- qbeta((1 - CI_level) / 2, alpha_n, beta_n)
  pred_upper <- qbeta((1 + CI_level) / 2, alpha_n, beta_n)

  # Output list
  prediction <- list(
    X        = X,
    Xnew     = Xnew,
    alpha_n  = alpha_n,
    beta_n   = beta_n,
    mean     = pred_mean,
    variance = pred_var,
    lower    = pred_lower,
    upper    = pred_upper,
    CI_level  = CI_level)

  # Posterior classification label (only for classification data)
  if (all(m == 1)) {
    prediction$class <- ifelse(pred_mean > threshold, 1, 0)
    prediction$threshold <- threshold
  }

  class(prediction) <- "predict_BKP"
  return(prediction)
}
