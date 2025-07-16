#' @name predict
#'
#' @title Predict method for Beta Kernel Process (BKP) models
#'
#' @description Generate predictions from a fitted BKP model at new input
#'   locations.
#'
#' @param object A fitted BKP model object returned by \code{\link{fit.BKP}}.
#' @param Xnew A matrix (or vector) of new input points at which to predict.
#' @param CI_size Confidence interval level (default = 0.05 for 95% CI).
#' @param threshold Classification threshold for posterior mean (default = 0.5).
#' @param ... Additional arguments passed to generic predict functions
#'   (currently not used, included for S3 method consistency).
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{X}{Original new input locations}
#'   \item{mean}{Posterior mean of success probability}
#'   \item{variance}{Posterior variance}
#'   \item{lower}{Lower bound of CI for success probability}
#'   \item{upper}{Upper bound of CI for success probability}
#'   \item{class}{(Optional) Predicted binary label if all trials are \( m_i = 1 \)}
#' }
#'
#' @author Jiangyan Zhao, Kunhai Qing, Jin Xu
#'
#' @keywords BKP
#'
#' @examples
#' ### 1D
#' set.seed(123)
#' n <- 30
#' Xbounds <- matrix(c(-2,2), nrow=1)
#' x <- tgp::lhs(n = n, rect = Xbounds)
#' m <- sample(100, n, replace = TRUE)
#' true_pi <- (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' y <- rbinom(n, size = m, prob = true_pi)
#' model1 <- fit.BKP(x, y, m, Xbounds=Xbounds)
#' xx = matrix(seq(-2, 2, length = 100), ncol=1) #new data points
#' head(predict(model1,xx))
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
#' m <- sample(100, n, replace = TRUE)
#' true_pi <- pnorm(f(x))
#' y <- rbinom(n, size = m, prob = true_pi)
#' model2 <- fit.BKP(x, y, m)
#' xx1 <- seq(Xbounds[1,1], Xbounds[1,2], length.out = 100)
#' xx2 <- seq(Xbounds[2,1], Xbounds[2,2], length.out = 100)
#' xx <- expand.grid(xx1 = xx1, xx2 = xx2)
#' head(predict(model2,xx))
#'
#' @export
#' @method predict BKP

predict.BKP <- function(object, Xnew, CI_size = 0.05, threshold = 0.5, ...)
{
  if (!inherits(object, "BKP")) {
    stop("The input is not of class 'BKP'. Please provide a model fitted with 'fit.BKP()'.")
  }

  BKPmodel <- object

  # Extract components
  Xnorm   <- BKPmodel$Xnorm
  y       <- BKPmodel$y
  m       <- BKPmodel$m
  theta   <- BKPmodel$theta_opt
  kernel  <- BKPmodel$kernel
  prior   <- BKPmodel$prior
  r0      <- BKPmodel$r0
  p0      <- BKPmodel$p0
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
  Xnew_norm <- sweep(Xnew, 2, Xbounds[, 1], "-")
  Xnew_norm <- sweep(Xnew_norm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

  # Compute kernel matrix
  K <- kernel_matrix(Xnew_norm, Xnorm, theta = theta, kernel = kernel) # m*n matrix

  # get the prior parameters: alpha0(x) and beta0(x)
  prior_par <- get_prior(prior = prior, r0 = r0, p0 = p0, y = y, m = m, K = K)
  alpha0 <- prior_par$alpha0
  beta0 <- prior_par$beta0

  # Posterior parameters
  alpha_n <- alpha0 + as.vector(K %*% y)
  beta_n  <- beta0 + as.vector(K %*% (m - y))

  # Predictive mean and variance
  pi_mean <- alpha_n / (alpha_n + beta_n)
  pi_var  <- pi_mean * (1 - pi_mean) / (alpha_n + beta_n + 1)

  # Confidence intervals
  pi_lower <- qbeta(CI_size / 2, alpha_n, beta_n)
  pi_upper <- qbeta(1 - CI_size / 2, alpha_n, beta_n)


  # Output
  prediction <- data.frame(Xnew,
                           mean = pi_mean,
                           variance = pi_var,
                           lower = pi_lower,
                           upper = pi_upper)
  names(prediction)[1:d] <- paste0("x", 1:d)

  # Posterior classification label (only for classification data)
  if (all(m == 1)) {
    prediction$class <- ifelse(pi_mean > threshold, 1, 0)
  }

  return(prediction)
}
