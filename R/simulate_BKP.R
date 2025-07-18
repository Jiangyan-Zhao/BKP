#' @name simulate
#'
#' @title Simulate from a Fitted BKP or DKP Model
#'
#' @description Generates random samples from the posterior distributions of a
#'   fitted BKP or DKP model at new input locations. For BKP models, optional
#'   binary classification labels can be returned based on a user-specified
#'   threshold. For DKP models, optional multiclass labels can be assigned using
#'   the maximum a posteriori (MAP) rule.
#'
#' @param object An object of class \code{"BKP"} or \code{"DKP"}, typically
#'   returned by \code{\link{fit.BKP}} or \code{\link{fit.DKP}}.
#' @param Xnew A numeric matrix or vector of new input locations for simulation.
#' @param n_sim Number of posterior samples to generate (default = \code{1}).
#' @param threshold Classification threshold for binary output (only used for
#'   BKP). If specified, the output will include binary class labels with values
#'   above the threshold classified as 1 (default is \code{NULL}).
#' @param seed Optional integer seed for reproducibility.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{\code{sims}}{
#'     BKP: A numeric matrix of dimension \code{nrow(Xnew) × n_sim} with simulated success probabilities.\cr
#'     DKP: A numeric array of dimension \code{n_sim × q × nrow(Xnew)} representing simulated class probabilities
#'     from Dirichlet distributions, where \code{q} is the number of classes.
#'   }
#'   \item{\code{mean}}{
#'     BKP: A numeric vector of posterior means of success probabilities.\cr
#'     DKP: A numeric matrix of posterior mean class probabilities (\code{nrow(Xnew) × q}).
#'   }
#'   \item{\code{class}}{
#'     BKP: A binary classification matrix (\code{nrow(Xnew) × n_sim}) if \code{threshold} is provided.\cr
#'     DKP: A vector of MAP class labels if training data is one-hot encoded (i.e., each row of \code{Y} sums to 1).
#'   }
#' }
#'
#' @seealso \code{\link{fit.BKP}}, \code{\link{fit.DKP}},
#'   \code{\link{predict.BKP}}, \code{\link{predict.DKP}}
#'
#' @keywords BKP
#'
#' @examples
#' ## -------------------- BKP Simulation Example --------------------
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
#' model <- fit.BKP(X, y, m, Xbounds=Xbounds)
#'
#' # Simulate 5 posterior draws of success probabilities
#' Xnew <- matrix(seq(-2, 2, length.out = 10), ncol = 1)
#' simulate(model, Xnew, n_sim = 5)
#'
#' # Simulate binary classifications (threshold = 0.5)
#' simulate(model, Xnew, n_sim = 5, threshold = 0.5)
#'
#' @export
#' @method simulate BKP

simulate.BKP <- function(object, Xnew, n_sim = 1, threshold = NULL, seed = NULL, ...) {
  if (!inherits(object, "BKP")) {
    stop("The input must be of class 'BKP'. Please provide a model fitted with 'fit.BKP()'.")
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

  # --- Check and normalize Xnew ---
  if (is.null(nrow(Xnew))) Xnew <- matrix(Xnew, nrow = 1)
  Xnew <- as.matrix(Xnew)
  if (ncol(Xnew) != d) {
    stop("Xnew must have the same number of columns as the training input.")
  }
  n_new <- nrow(Xnew)
  Xnew_norm <- sweep(Xnew, 2, Xbounds[, 1], "-")
  Xnew_norm <- sweep(Xnew_norm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

  # --- Compute kernel matrix between Xnew and training X ---
  K <- kernel_matrix(Xnew_norm, Xnorm, theta = theta, kernel = kernel)

  # --- Get prior parameters ---
  prior_par <- get_prior(prior = prior, r0 = r0, p0 = p0, y = y, m = m, K = K)
  alpha0 <- prior_par$alpha0
  beta0 <- prior_par$beta0

  # --- Compute posterior Beta parameters ---
  alpha_n <- pmax(alpha0 + as.vector(K %*% y), 1e-10)
  beta_n  <- pmax(beta0 + as.vector(K %*% (m - y)), 1e-10)

  # --- Simulate from posterior Beta distributions ---
  if (!is.null(seed)) set.seed(seed)
  sims <- matrix(rbeta(n_new * n_sim,
                       shape1 = rep(alpha_n, n_sim),
                       shape2 = rep(beta_n, n_sim)),
                 nrow = n_new, ncol = n_sim)
  colnames(sims) <- paste0("sim", 1:n_sim)
  rownames(sims) <- paste0("x", 1:n_new)

  # --- Optional: binary classification ---
  class_pred <- NULL
  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0 || threshold >= 1) {
      stop("Threshold must be a numeric value strictly between 0 and 1.")
    }
    class_pred <- ifelse(sims > threshold, 1L, 0L)
  }

  # --- Posterior mean ---
  pi_mean <- alpha_n / (alpha_n + beta_n)

  return(list(
    sims = sims,        # [n_new × n_sim]: simulated probabilities
    mean = pi_mean,     # [n_new]: posterior mean
    class = class_pred  # [n_new × n_sim]: binary labels (if threshold provided)
  ))
}


#' @export
simulate <- function(object, ...) {
  UseMethod("simulate")
}
