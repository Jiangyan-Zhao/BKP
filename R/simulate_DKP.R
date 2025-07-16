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
#' true_pi <- matrix(c(true_pi/2,true_pi/2,1-true_pi),nrow = n, byrow = FALSE)
#' m <- sample(100, n, replace = TRUE)
#' Y <- matrix(0, nrow = n, ncol = 3)
#' for (i in 1:n) {
#'   Y[i, ] <- rmultinom(n=1, size=m[i], prob=true_pi[i, ])
#' }
#' DKPmodel <- fit.DKP(x, Y, Xbounds = Xbounds)
#' simulate(DKPmodel,Xnew = 0.5, n_sim = 5)
#'
#' @export
#' @method simulate DKP

simulate.DKP <- function(object, Xnew, n_sim = 1, seed = NULL, ...) {
  if (!inherits(object, "DKP")) {
    stop("The input must be of class 'DKP'. Please provide a model fitted with 'fit.DKP()'.")
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
  q       <- ncol(Y)

  # Ensure Xnew is a matrix
  if (is.null(nrow(Xnew))) Xnew <- matrix(Xnew, nrow = 1)
  Xnew <- as.matrix(Xnew)
  if (ncol(Xnew) != d) {
    stop("Xnew must have the same number of columns as the training input.")
  }
  n_new <- nrow(Xnew)

  # Normalize Xnew to [0,1]^d
  Xnew_norm <- sweep(Xnew, 2, Xbounds[, 1], "-")
  Xnew_norm <- sweep(Xnew_norm, 2, Xbounds[, 2] - Xbounds[, 1], "/")


  # Kernel matrix
  K <- kernel_matrix(Xnew_norm, Xnorm, theta = theta, kernel = kernel)

  # Prior parameters
  alpha0 <- get_prior_dkp(prior = prior, r0 = r0, p0 = p0, Y = Y, K = K)

  # Compute posterior alpha
  alpha_n <- as.matrix(alpha0) + as.matrix(K %*% Y)
  alpha_n <- pmax(alpha_n, 1e-10) # Avoid numerical issues

  # Simulate posterior samples: array [n_sim × q × n_new]
  if (!is.null(seed)) set.seed(seed)
  sims <- array(0, dim = c(n_sim, q, n_new))
  for (i in 1:n_new) {
    shape_mat <- matrix(rgamma(n_sim * q, shape = rep(alpha_n[i, ], each = n_sim), rate = 1),
                        nrow = n_sim, byrow = FALSE)
    sims[,,i] <- shape_mat / rowSums(shape_mat)
  }
  # Name dimensions
  dimnames(sims) <- list(
    paste0("sim", 1:n_sim),
    paste0("Class", 1:q),
    paste0("Xnew_", 1:n_new)
  )

  # Compute posterior mean
  pi_mean <- alpha_n / rowSums(alpha_n)

  # Optional: hard class prediction
  class_pred <- NULL
  if (all(rowSums(Y) == 1)) {
    class_pred <- max.col(pi_mean)
  }

  # Return results
  return(list(
    sims = sims,             # [n_sim × q × n_new]
    mean = pi_mean,          # [n_new × q]
    class = class_pred       # [n_new] — no comma here!
  ))
}

