#' @name simulate
#'
#' @title Simulate method for Beta Kernel Process (BKP) models
#'
#' @description Generate random samples from the posterior Beta distributions
#'   of a fitted BKP model at new input locations. Optionally return
#'   classification labels based on a threshold.
#'
#' @param object A fitted BKP model object from \code{\link{fit.BKP}}.
#' @param Xnew A matrix (or vector) of new input points at which to simulate.
#' @param n_sim Number of simulation replicates (default = 1).
#' @param threshold Classification threshold for labeling (default = NULL;
#'   if provided, returns 0/1 labels instead of probabilities).
#' @param seed Optional random seed for reproducibility.
#' @param ... Additional arguments (currently unused).
#'
#' @return A matrix of dimension \code{nrow(Xnew) x n_sim}:
#'   \itemize{
#'     \item If \code{threshold = NULL}, each column contains simulated success probabilities.
#'     \item If \code{threshold} is provided, the matrix contains binary class labels (0/1).
#'   }
#'
#' @author Jiangyan Zhao, Kunhai Qing, Jin Xu
#'
#' @keywords BKP
#'
#' @examples
#' ## 1D example: simulate probabilities and classifications
#' set.seed(123)
#' n <- 30
#' x <- seq(-2, 2, length = n)
#' true_pi <- (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' m <- sample(50:100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#' df <- data.frame(x = x, y = y, m = m)
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' model <- fit.BKP(df, Xbounds = Xbounds)
#'
#' # Simulate 5 posterior draws of success probabilities at new points
#' Xnew <- matrix(seq(-2, 2, length.out = 10), ncol = 1)
#' simulate(model, Xnew, n_sim = 5)
#'
#' # Simulate 5 posterior binary classifications using threshold 0.5
#' simulate(model, Xnew, n_sim = 5, threshold = 0.5)
#'
#' @export
#' @method simulate BKP

simulate.BKP <- function(object, Xnew, n_sim = 1, threshold = NULL, seed = NULL, ...) {
  if (!inherits(object, "BKP")) {
    stop("The input must be of class 'BKP'. Please provide a model fitted with 'fit.BKP()'.")
  }

  if (!is.null(seed)) set.seed(seed)

  BKPmodel <- object

  # Extract components
  Xnorm   <- BKPmodel$Xnorm
  y       <- BKPmodel$y
  m       <- BKPmodel$m
  theta   <- BKPmodel$bestTheta
  kernel  <- BKPmodel$kernel
  prior   <- BKPmodel$prior
  r0      <- BKPmodel$r0
  p0      <- BKPmodel$p0
  Xbounds <- BKPmodel$Xbounds
  d       <- ncol(Xnorm)

  # Ensure Xnew is a matrix
  if (is.null(nrow(Xnew))) Xnew <- matrix(Xnew, nrow = 1)
  Xnew <- as.matrix(Xnew)
  if (ncol(Xnew) != d) {
    stop("Xnew must have the same number of columns as the training input.")
  }

  # Normalize Xnew to [0,1]^d
  Xnew_norm <- sweep(Xnew, 2, Xbounds[, 1], "-")
  Xnew_norm <- sweep(Xnew_norm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

  # Kernel matrix
  K <- kernel_matrix(Xnew_norm, Xnorm, theta = theta, kernel = kernel)

  # Prior parameters
  prior_par <- get_prior(prior = prior, r0 = r0, p0 = p0, y = y, m = m, K = K)
  alpha0 <- prior_par$alpha0
  beta0 <- prior_par$beta0

  # Posterior parameters
  alpha_n <- alpha0 + as.vector(K %*% y)
  beta_n  <- beta0 + as.vector(K %*% (m - y))

  # Simulate from Beta
  sims <- matrix(rbeta(length(alpha_n) * n_sim,
                       shape1 = rep(alpha_n, n_sim),
                       shape2 = rep(beta_n, n_sim)),
                 nrow = length(alpha_n), ncol = n_sim)

  # Optional: convert to binary class labels
  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0 || threshold >= 1) {
      stop("Threshold must be a number between 0 and 1.")
    }
    sims <- ifelse(sims > threshold, 1L, 0L)
  }

  colnames(sims) <- paste0("sim", 1:n_sim)
  return(sims)
}


#' @export
simulate <- function(object, ...) {
  UseMethod("simulate")
}
