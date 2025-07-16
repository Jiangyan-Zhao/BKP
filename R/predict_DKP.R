#' @name predict
#'
#' @title Predict method for Dirichlet Kernel Process (DKP) models
#'
#' @description
#' Generate predictions from a fitted DKP model at new input locations.
#'
#' @param object A fitted DKP model object returned by \code{\link{fit.DKP}}.
#' @param Xnew A matrix (or vector) of new input points at which to predict.
#' @param CI_size Confidence interval level (default = 0.05 for 95% CI).
#' @param ... Additional arguments passed to generic predict functions (currently not used, included for S3 method consistency).
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
#' @keywords DKP
#'
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
#' predict(DKPmodel,Xnew = 0.5)
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
#' predict(DKPmodel, Xnew = c(0.5,0.5))
#'
#' @export
#' @method predict DKP
#' @importFrom stats qbeta


predict.DKP <- function(object, Xnew, CI_size = 0.05, ...)
{
  if (!inherits(object, "DKP")) {
    stop("The input is not of class 'DKP'. Please provide a model fitted with 'fit.DKP()'.")
  }

  DKPmodel <- object

  # Extract components
  Xnorm   <- DKPmodel$Xnorm
  Y       <- DKPmodel$Y
  theta   <- DKPmodel$theta_opt
  kernel  <- DKPmodel$kernel
  prior   <- DKPmodel$prior
  r0      <- DKPmodel$r0
  p0      <- DKPmodel$p0
  Xbounds <- DKPmodel$Xbounds
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
  alpha0 <- get_prior_dkp(prior = prior, r0 = r0, p0 = p0, Y = Y, K = K)

  # Posterior parameters
  alpha_n <- alpha0 + as.matrix(K %*% Y)

  # Predictive mean and variance
  pi_mean <- alpha_n / rowSums(alpha_n)
  pi_var  <- alpha_n * (rowSums(alpha_n) - alpha_n) / (rowSums(alpha_n)^2 * (rowSums(alpha_n) + 1))

  # Confidence intervals
  pi_lower <- qbeta(CI_size / 2, alpha_n, rowSums(alpha_n) - alpha_n)
  pi_upper <- qbeta(1 - CI_size / 2, alpha_n, rowSums(alpha_n) - alpha_n)

  # Output
  prediction <- data.frame(Xnew,
                           mean = pi_mean,
                           variance = pi_var,
                           lower = pi_lower,
                           upper = pi_upper)
  names(prediction)[1:d] <- paste0("x", 1:d)

  # Posterior classification label (only for classification data)
  if (all(rowSums(Y) == 1)) {
    prediction$class <-max.col(pi_mean)
  }

  return(prediction)
}
