#' @name simulate
#'
#' @title Simulate from Fitted BKP Package Models
#'
#' @description Generate random posterior draws from fitted \code{BKP},
#'   \code{DKP}, \code{TwinBKP}, and \code{TwinDKP} model objects at training
#'   or new input locations.
#'
#'   For \code{BKP} and \code{TwinBKP} models, posterior samples are generated
#'   from Beta distributions for latent success probabilities. Optionally,
#'   binary class labels can be derived by applying a user-specified
#'   classification threshold.
#'
#'   For \code{DKP} and \code{TwinDKP} models, posterior samples are generated
#'   from Dirichlet distributions for latent class-probability vectors. If the
#'   training responses are single-label, that is, one-hot encoded, class labels
#'   are additionally assigned using the maximum posterior simulated probability.
#'
#' @param object A fitted model object of class \code{"BKP"}, \code{"DKP"},
#'   \code{"TwinBKP"}, or \code{"TwinDKP"}, typically returned by
#'   \code{\link{fit_BKP}}, \code{\link{fit_DKP}},
#'   \code{\link{fit_TwinBKP}}, or \code{\link{fit_TwinDKP}}.
#' @param Xnew A numeric vector or matrix of new input locations at which
#'   posterior simulations are generated. If \code{NULL}, simulations are
#'   generated at the training input locations.
#' @param nsim Number of posterior samples to generate. The default is
#'   \code{1}.
#' @param threshold Classification threshold for binary decisions, used for
#'   \code{BKP} and \code{TwinBKP} models. When specified, posterior draws
#'   exceeding the threshold are classified as 1, and those below or equal to
#'   the threshold are classified as 0. The default is \code{NULL}.
#' @param seed Optional integer seed for reproducibility.
#' @param ... Additional arguments passed to the corresponding \code{predict}
#'   method when \code{Xnew} is supplied.
#'
#' @return A list containing posterior simulations:
#' \describe{
#'   \item{\code{samples}}{
#'     For \code{BKP} and \code{TwinBKP}, a numeric matrix with one row per
#'     simulation location and one column per posterior draw. Each entry is a
#'     simulated latent success probability.\cr
#'     For \code{DKP} and \code{TwinDKP}, a numeric array with dimensions
#'     simulation locations \eqn{\times} classes \eqn{\times} posterior draws.
#'     Each slice contains simulated latent class probabilities from the
#'     posterior Dirichlet distribution.
#'   }
#'   \item{\code{mean}}{
#'     Posterior mean at the simulation locations. For \code{BKP} and
#'     \code{TwinBKP}, this is a numeric vector of success-probability means.
#'     For \code{DKP} and \code{TwinDKP}, this is a matrix of class-probability
#'     means.
#'   }
#'   \item{\code{class}}{
#'     Simulated class labels when applicable. For \code{BKP} and
#'     \code{TwinBKP}, this is returned when \code{threshold} is supplied. For
#'     \code{DKP} and \code{TwinDKP}, this is returned when the training
#'     response is one-hot encoded.
#'   }
#'   \item{\code{X}}{The original training input locations.}
#'   \item{\code{Xnew}}{The new input locations used for simulation. If
#'     \code{NULL}, simulations are returned at the training input locations.}
#'   \item{\code{threshold}}{The binary classification threshold, returned for
#'     \code{BKP} and \code{TwinBKP} simulations when supplied.}
#'   \item{\code{ess}}{Effective-sample-size calibration method inherited from
#'     the fitted model. For \code{TwinBKP} and \code{TwinDKP}, this is
#'     \code{"none"}.}
#'   \item{\code{ess_info}}{ESS diagnostics for \code{BKP} and \code{DKP}
#'     simulations when available. \code{TwinBKP} and \code{TwinDKP} use
#'     uncalibrated global-local posterior updates and return
#'     \code{ess_info = NULL}.}
#' }
#'
#' @seealso \code{\link{fit_BKP}}, \code{\link{fit_DKP}},
#'   \code{\link{fit_TwinBKP}}, and \code{\link{fit_TwinDKP}} for model fitting;
#'   \code{\link{predict.BKP}}, \code{\link{predict.DKP}},
#'   \code{\link{predict.TwinBKP}}, and \code{\link{predict.TwinDKP}} for
#'   posterior prediction.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}. arXiv. <doi:10.48550/arXiv.2508.10447>.
#'
#' @keywords BKP
#'
#' @examples
#' # -------------------------- BKP and TwinBKP ---------------------------
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
#' model <- fit_BKP(X, y, m, Xbounds=Xbounds)
#'
#' # Simulate 5 posterior draws of success probabilities
#' Xnew <- matrix(seq(-2, 2, length.out = 5), ncol = 1)
#' simulate(model, Xnew = Xnew, nsim = 5)
#'
#' # Simulate binary classifications (threshold = 0.5)
#' simulate(model, Xnew = Xnew, nsim = 5, threshold = 0.5)
#'
#'#' \dontrun{
#' # Larger TwinBKP example
#' n <- 1000
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit TwinBKP model
#' model <- fit_TwinBKP(X, y, m, Xbounds=Xbounds)
#'
#' # Simulate 5 posterior draws of success probabilities
#' simulate(model, Xnew = Xnew, nsim = 5)
#'
#' # Simulate binary classifications (threshold = 0.5)
#' simulate(model, Xnew = Xnew, nsim = 5, threshold = 0.5)
#' }
#'
#' @export
#' @method simulate BKP

simulate.BKP <- function(object, nsim = 1, seed = NULL, Xnew = NULL, threshold = NULL, ...)
{
  # ---------------- Argument Checking ----------------
  if (!is.numeric(nsim) || length(nsim) != 1L ||
      is.na(nsim) || !is.finite(nsim) ||
      nsim <= 0 || nsim != as.integer(nsim)) {
    stop("`nsim` must be a positive integer.")
  }
  nsim <- as.integer(nsim)

  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L ||
       is.na(seed) || !is.finite(seed) ||
       seed != as.integer(seed))) {
    stop("`seed` must be a single integer or NULL.")
  }

  d <- ncol(object$X)

  if (!is.null(Xnew)) {
    if (is.null(dim(Xnew))) {
      if (d == 1L) {
        Xnew <- matrix(Xnew, ncol = 1L)
      } else {
        Xnew <- matrix(Xnew, nrow = 1L)
      }
    } else {
      Xnew <- as.matrix(Xnew)
    }

    if (!is.numeric(Xnew)) {
      stop("'Xnew' must be numeric.")
    }
    if (nrow(Xnew) < 1L || ncol(Xnew) < 1L) {
      stop("'Xnew' must have at least one row and one column.")
    }
    if (ncol(Xnew) != d) {
      stop("The number of columns in 'Xnew' must match the original input dimension.")
    }
    if (anyNA(Xnew) || any(!is.finite(Xnew))) {
      stop("'Xnew' must contain only finite values with no NA, NaN, or Inf.")
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
    prediction <- predict.BKP(object, Xnew = Xnew, type = "probability", ...)
    alpha_n <- prediction$alpha_n
    beta_n <- prediction$beta_n
    ess_info <- prediction$ess_info
  } else {
    # Use training data
    alpha_n <- object$alpha_n
    beta_n  <- object$beta_n
    ess_info <- object$ess_info
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
    threshold = threshold,  # classification threshold
    ess       = if (is.null(object$ess)) "none" else object$ess,
    ess_info  = ess_info
  )

  class(simulation) <- "simulate_BKP"
  return(simulation)
}
