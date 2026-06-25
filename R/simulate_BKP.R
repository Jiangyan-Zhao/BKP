#' @name simulate
#'
#' @title Simulate from a Fitted BKP, DKP, or TwinBKP Model
#'
#' @description Generates random draws from the posterior distribution of
#'   fitted BKP, DKP, or TwinBKP models at training or new input locations.
#'
#'   For BKP and TwinBKP models, posterior samples are generated from Beta
#'   distributions characterizing success probabilities. Optionally, binary
#'   class labels can be derived by applying a user-specified classification
#'   threshold.
#'
#'   For DKP models, posterior samples are generated from Dirichlet
#'   distributions characterizing class probabilities. If training responses are
#'   single-label (i.e., one-hot encoded), class labels may additionally be
#'   assigned using the maximum a posteriori (MAP) rule.
#'
#' @param object An object of class \code{"BKP"}, \code{"DKP"}, or
#'   \code{"TwinBKP"}, typically returned by \code{\link{fit_BKP}},
#'   \code{\link{fit_DKP}}, or \code{\link{fit_TwinBKP}}.
#' @param Xnew A numeric matrix or vector of new input locations at which
#'   simulations are generated.
#' @param nsim Number of posterior samples to generate (default = \code{1}).
#' @param threshold Classification threshold for binary decisions, used for
#'   BKP and TwinBKP models. When specified, posterior draws exceeding the
#'   threshold are classified as 1, and those below as 0. The default is
#'   \code{NULL}.
#' @param seed Optional integer seed for reproducibility.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{samples}}{
#'     For \strong{BKP} and \strong{TwinBKP}: A numeric matrix of size
#'     \code{nrow(Xnew) × nsim}, where each column corresponds to one posterior
#'     draw of success probabilities.\cr
#'     For \strong{DKP}: A numeric array of dimension
#'     \code{nrow(Xnew) × q × nsim}, containing simulated class probabilities
#'     from Dirichlet posteriors, where \code{q} is the number of classes.
#'   }
#'
#'   \item{\code{mean}}{
#'     For \strong{BKP} and \strong{TwinBKP}: A numeric vector of posterior mean
#'     success probabilities at each input location.\cr
#'     For \strong{DKP}: A numeric matrix of dimension \code{nrow(Xnew) × q},
#'     containing posterior mean class probabilities.
#'   }
#'
#'   \item{\code{class}}{
#'     For \strong{BKP} and \strong{TwinBKP}: An integer matrix of dimension
#'     \code{nrow(Xnew) × nsim}, indicating simulated binary class labels
#'     \code{0}/\code{1}, returned when \code{threshold} is specified.\cr
#'     For \strong{DKP}: An integer matrix of dimension
#'     \code{nrow(Xnew) × nsim}, where each entry corresponds to a MAP-predicted
#'     class label, returned only when training data is single-label.
#'   }
#'
#'   \item{\code{X}}{The training input matrix used to fit the BKP, DKP, or
#'     TwinBKP model.}
#'
#'   \item{\code{Xnew}}{The new input locations at which simulations are generated.}
#'
#'   \item{\code{threshold}}{The classification threshold used for generating
#'     binary class labels (if provided).}
#'
#'   \item{\code{ess}}{Effective-sample-size calibration method inherited from
#'     the fitted model, returned for BKP objects when available.}
#'
#'   \item{\code{ess_info}}{ESS diagnostics for BKP simulations when available.
#'     TwinBKP uses the uncalibrated global-local posterior update and does not
#'     return ESS diagnostics.}
#' }
#'
#' @seealso \code{\link{fit_BKP}}, \code{\link{fit_DKP}}, and
#'   \code{\link{fit_TwinBKP}} for model fitting; \code{\link{predict.BKP}},
#'   \code{\link{predict.DKP}}, and \code{\link{predict.TwinBKP}} for posterior
#'   prediction.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}. arXiv. \doi{10.48550/arXiv.2508.10447}
#'
#' @keywords BKP DKP TwinBKP
#'
#' @examples
#' ## -------------------- BKP --------------------
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
#' @export
#' @method simulate BKP

simulate.BKP <- function(object, nsim = 1, seed = NULL, Xnew = NULL, threshold = NULL, ...)
{
  # ---------------- Argument Checking ----------------
  if (!is.numeric(nsim) || length(nsim) != 1 || nsim <= 0 || nsim != as.integer(nsim)) {
    stop("`nsim` must be a positive integer.")
  }
  nsim <- as.integer(nsim)

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1 || seed != as.integer(seed))) {
    stop("`seed` must be a single integer or NULL.")
  }

  d <- ncol(object$X)
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
