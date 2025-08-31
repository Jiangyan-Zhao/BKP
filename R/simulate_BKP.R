#' @name simulate
#'
#' @title Simulate from a Fitted BKP or DKP Model
#'
#' @description Generates random samples from the posterior predictive
#'   distribution of a fitted BKP or DKP model at new input locations.
#'
#'   For BKP models, posterior samples are drawn from Beta distributions
#'   representing success probabilities, with optional binary class labels
#'   determined by a threshold.
#'
#'   For DKP models, posterior samples are drawn from Dirichlet distributions
#'   representing class probabilities, with optional class labels determined by
#'   the maximum a posteriori (MAP) rule if training responses are one-hot
#'   encoded.
#'
#' @param object An object of class \code{"BKP"} or \code{"DKP"}, typically
#'   returned by \code{\link{fit.BKP}} or \code{\link{fit.DKP}}.
#' @param Xnew A numeric matrix or vector of new input locations for simulation.
#' @param nsim Number of posterior samples to generate (default = \code{1}).
#' @param threshold Classification threshold for binary output (only used for
#'   BKP). If specified, the output will include binary class labels with values
#'   above the threshold classified as 1 (default is \code{NULL}).
#' @param seed Optional integer seed for reproducibility.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{samples}}{
#'     For \strong{BKP}: A numeric matrix of dimension \code{nrow(Xnew) × nsim}, containing simulated success probabilities.\cr
#'     For \strong{DKP}: A numeric array of dimension \code{nsim × q × nrow(Xnew)}, containing simulated class probabilities
#'     from Dirichlet posteriors, where \code{q} is the number of classes.
#'   }
#'
#'   \item{\code{mean}}{
#'     For \strong{BKP}: A numeric vector of posterior mean success probabilities at each \code{Xnew}.\cr
#'     For \strong{DKP}: A numeric matrix of dimension \code{nrow(Xnew) × q}, containing posterior mean class probabilities.
#'   }
#'
#'   \item{\code{class}}{
#'     For \strong{BKP}: A binary matrix of dimension \code{nrow(Xnew) × nsim} indicating simulated class labels (0/1),
#'     returned if \code{threshold} is specified.\cr
#'     For \strong{DKP}: A numeric matrix of dimension \code{nrow(Xnew) × nsim} containing MAP-predicted class labels,
#'     returned only when training data is single-label (i.e., each row of \code{Y} sums to 1).
#'   }
#' }
#'
#' @seealso \code{\link{fit.BKP}}, \code{\link{fit.DKP}},
#'   \code{\link{predict.BKP}}, \code{\link{predict.DKP}}
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv.
#'   https://doi.org/10.48550/arXiv.2508.10447.
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
#' Xnew <- matrix(seq(-2, 2, length.out = 5), ncol = 1)
#' simulate(model, Xnew = Xnew, nsim = 5)
#'
#' # Simulate binary classifications (threshold = 0.5)
#' simulate(model, Xnew = Xnew, nsim = 5, threshold = 0.5)
#'
#' @export
#' @method simulate BKP

simulate.BKP <- function(object, nsim = 1, seed = NULL, ..., Xnew = NULL, threshold = NULL)
{
  # ---------------- Argument Checking ----------------
  if (!is.numeric(nsim) || length(nsim) != 1 || nsim <= 0 || nsim != as.integer(nsim)) {
    stop("`nsim` must be a positive integer.")
  }
  nsim <- as.integer(nsim)

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1 || seed != as.integer(seed))) {
    stop("`seed` must be a single integer or NULL.")
  }

  d <- ncol(object$Xnorm)
  if (!is.null(Xnew)) {
    if (is.null(nrow(Xnew))) Xnew <- matrix(Xnew, nrow = 1)
    if (!is.numeric(Xnew)) {
      stop("`Xnew` must be numeric.")
    }
    Xnew <- as.matrix(Xnew)
    if (ncol(Xnew) != d) {
      stop(sprintf("`Xnew` must have %d columns, matching the training inputs.", d))
    }
  }

  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0 || threshold >= 1) {
      stop("`threshold` must be a numeric value strictly between 0 and 1 (e.g., 0.5).")
    }
  }

  # ---------------- Core Computation ----------------
  if (!is.null(seed)) set.seed(seed)

  if (!is.null(Xnew)) {
    # complete posterior parameters at new inputs
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

    # --- Normalize new inputs ---
    Xnew_norm <- sweep(Xnew, 2, Xbounds[, 1], "-")
    Xnew_norm <- sweep(Xnew_norm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

    # --- Compute kernel matrix between Xnew and training X ---
    K <- kernel_matrix(Xnew_norm, Xnorm, theta = theta, kernel = kernel)

    # --- Get prior parameters ---
    prior_par <- get_prior(prior = prior, model = "BKP",
                           r0 = r0, p0 = p0, y = y, m = m, K = K)
    alpha0 <- prior_par$alpha0
    beta0 <- prior_par$beta0

    # --- Compute posterior Beta parameters ---
    alpha_n <- pmax(alpha0 + as.vector(K %*% y), 1e-10)
    beta_n  <- pmax(beta0 + as.vector(K %*% (m - y)), 1e-10)
  }else{
    # Use training data
    alpha_n <- pmax(object$alpha_n, 1e-10)
    beta_n  <- pmax(object$beta_n, 1e-10)
  }

  # --- Simulate from posterior Beta distributions ---
  n_new <- ifelse(!is.null(Xnew), nrow(Xnew), nrow(object$X))
  samples <- matrix(rbeta(n_new * nsim,
                       shape1 = rep(alpha_n, nsim),
                       shape2 = rep(beta_n, nsim)),
                 nrow = n_new, ncol = nsim)
  colnames(samples) <- paste0("sim", 1:nsim)
  rownames(samples) <- paste0("x", 1:n_new)

  # --- Optional: binary classification ---
  class_pred <- NULL
  if (!is.null(threshold)) {
    class_pred <- ifelse(samples > threshold, 1L, 0L)
  }

  # --- Posterior mean ---
  pi_mean <- alpha_n / (alpha_n + beta_n)

  simulation <- list(
    samples   = samples,    # [n_new × nsim]: simulated probabilities
    mean      = pi_mean,    # [n_new]: posterior mean
    class     = class_pred, # [n_new × nsim]: binary labels (if threshold provided)
    X         = object$X,   # [n × d]: training inputs
    Xnew      = Xnew,       # [n_new × d]: new inputs (if provided)
    threshold = threshold   # classification threshold
  )

  class(simulation) <- "simulate.BKP"
  return(simulation)
}
