#' @name simulate
#'
#' @title Simulate method for Dirichlet Kernel Process (DKP) models
#'
#' @description Generate random samples from the posterior Dirichlet distributions
#'   of a fitted DKP model at new input locations. Optionally return
#'   classification labels based on a threshold.
#'
#' @param object A fitted DKP model object from \code{\link{fit.DKP}}.
#' @param Xnew A matrix (or vector) of new input points at which to simulate.
#' @param n_sim Number of simulation replicates (default = 1).
#' @param threshold Classification threshold for labeling (default = NULL)
#' @param seed Optional random seed for reproducibility.
#' @param ... Additional arguments (currently unused).
#'
#' @return A array of dimension \code{n_sim x q x nrow(Xnew)}:
#'   \itemize{
#'     \item If \code{threshold = NULL}, for every fold, each row contains simulated success probabilities.
#'     \item If \code{threshold} is provided, every fold contains multinomial class labels .
#'   }
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
#' simulate(DKPmodel,Xnew = 0.5, n_sim = 5)
#'
#' @export
#' @method simulate DKP

simulate.DKP <- function(object, Xnew, n_sim = 1, threshold = NULL, seed = NULL, ...) {
  if (!inherits(object, "DKP")) {
    stop("The input must be of class 'DKP'. Please provide a model fitted with 'fit.DKP()'.")
  }

  if (!is.null(seed)) set.seed(seed)

  DKPmodel <- object

  # Extract components
  Xnorm   <- DKPmodel$Xnorm
  Y       <- DKPmodel$Y
  theta   <- DKPmodel$bestTheta
  kernel  <- DKPmodel$kernel
  prior   <- DKPmodel$prior
  r0      <- DKPmodel$r0
  p0      <- DKPmodel$p0
  Xbounds <- DKPmodel$Xbounds
  d       <- ncol(Xnorm)
  q       <- ncol(Y)

  # Ensure Xnew is a matrix
  if (is.null(nrow(Xnew))) Xnew <- matrix(Xnew, nrow = 1)
  Xnew <- as.matrix(Xnew)
  if (ncol(Xnew) != d) {
    stop("Xnew must have the same number of columns as the training input.")
  }

  # Normalize Xnew to [0,1]^d
  Xnew_norm <- sweep(Xnew, 2, Xbounds[, 1], "-")
  Xnew_norm <- sweep(Xnew_norm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

  n_new       <- nrow(Xnew)
  # Kernel matrix
  K <- kernel_matrix(Xnew_norm, Xnorm, theta = theta, kernel = kernel)

  # Prior parameters
  prior_par <- get_prior_dkp(prior = prior, r0 = r0, p0 = p0, Y = Y, K = K)
  alpha0 <- prior_par$alpha0
  if (prior == "noninformative" || prior == "fixed"){
    alpha0 <- matrix(rep(alpha0,n_new),nrow = n_new, byrow = T)
  }
  # Compute posterior alpha
  alpha_n <- as.matrix(alpha0) + as.matrix(K %*% Y)

  # Numerical stabilization: avoid log(0) or NaNs
  alpha_n <- pmax(alpha_n, 1e-10)

  # rdirichlet
  rdirichlet <- function(n, alpha) {
    k <- length(alpha)
    x <- matrix(rgamma(n * k, shape = alpha, rate = 1), nrow = n, byrow = TRUE)
    samples <- x / rowSums(x)
    return(samples)
  }
  # Simulate from Beta
  sims <- array(0,dim = c(n_sim,q,n_new))
  for (i in 1:n_new ) {
    sims[,,i] <- rdirichlet(n_sim, alpha_n[i,])
  }

  # Optional: convert to multinomial class labels
  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0 || threshold >= 1) {
      stop("Threshold must be a number between 0 and 1.")
    }
    sims <- ifelse(sims > threshold, 1L, 0L)
  }
  dimnames(sims) <- list(paste0("sim", 1:n_sim),
                         paste0("Ylabel", 1:q),
                         paste0("Xnew", 1:n_new))
  return(sims)
}

