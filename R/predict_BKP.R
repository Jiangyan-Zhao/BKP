#' @name predict
#'
#' @title Posterior Prediction for BKP Package Model Objects
#'
#' @description Generate posterior summaries from fitted \code{BKP},
#'   \code{DKP}, \code{TwinBKP}, and \code{TwinDKP} model objects at training or
#'   new input locations. For BKP and TwinBKP models, summaries can be returned
#'   either for the latent success probability or for a future success count
#'   under the Beta-Binomial posterior predictive distribution. For DKP and
#'   TwinDKP models, summaries can be returned either for the latent class
#'   probability vector or for future class counts using marginal Beta-Binomial
#'   posterior predictive distributions.
#'
#' @param object A fitted model object of class \code{"BKP"}, \code{"DKP"},
#'   \code{"TwinBKP"}, or \code{"TwinDKP"}, typically returned by
#'   \code{\link{fit_BKP}}, \code{\link{fit_DKP}},
#'   \code{\link{fit_TwinBKP}}, or \code{\link{fit_TwinDKP}}.
#' @param Xnew A numeric vector or matrix of new input locations at which to
#'   generate predictions. If \code{NULL}, predictions are returned for the
#'   training data.
#' @param CI_level Numeric between 0 and 1 specifying the credible level for
#'   posterior intervals (default \code{0.95} for 95% credible interval).
#' @param threshold Numeric between 0 and 1 specifying the classification
#'   threshold for binary predictions based on posterior mean, used for
#'   BKP and TwinBKP models. Default is \code{0.5}.
#' @param type Character string specifying the prediction target. The default
#'   \code{"probability"} returns posterior summaries for latent success
#'   probabilities (BKP and TwinBKP) or latent class probabilities (DKP). The
#'   option \code{"count"} returns posterior predictive summaries for future
#'   counts under the Beta-Binomial distribution (BKP and TwinBKP) or marginal
#'   Beta-Binomial distributions for each multinomial class (DKP).
#' @param Mnew Positive integer trial size used when \code{type = "count"}. It
#'   can be either a scalar, applied to all prediction points, or a vector with
#'   the same length as the number of prediction points. If \code{Xnew = NULL}
#'   and \code{Mnew = NULL}, the training trial sizes \code{object$m} (BKP and
#'   TwinBKP) or \code{rowSums(object$Y)} (DKP) are used.
#' @param ... Additional arguments passed to generic \code{predict} methods
#'   (currently not used; included for S3 method consistency).
#'
#' @return A list containing posterior or posterior predictive summaries:
#' \describe{
#'   \item{\code{X}}{The original training input locations.}
#'   \item{\code{Xnew}}{The input locations at which predictions are returned.
#'     If \code{Xnew = NULL}, predictions are returned at the training input
#'     locations.}
#'   \item{\code{alpha_n}}{Posterior shape or concentration parameters:
#'     \itemize{
#'       \item For \code{BKP} and \code{TwinBKP}, a vector of posterior Beta
#'       \eqn{\alpha} shape parameters.
#'       \item For \code{DKP} and \code{TwinDKP}, a matrix of posterior
#'       Dirichlet concentration parameters, with rows corresponding to
#'       prediction locations and columns corresponding to classes.
#'     }}
#'   \item{\code{beta_n}}{Posterior Beta \eqn{\beta} shape parameters, returned
#'     for \code{BKP} and \code{TwinBKP} objects.}
#'   \item{\code{mean}}{Mean of the prediction target:
#'     \itemize{
#'       \item For \code{BKP} and \code{TwinBKP} with
#'       \code{type = "probability"}, the posterior mean of the latent success
#'       probability.
#'       \item For \code{BKP} and \code{TwinBKP} with \code{type = "count"},
#'       the posterior predictive mean of a future success count.
#'       \item For \code{DKP} and \code{TwinDKP} with
#'       \code{type = "probability"}, a matrix of posterior mean class
#'       probabilities.
#'       \item For \code{DKP} and \code{TwinDKP} with \code{type = "count"},
#'       a matrix of marginal posterior predictive mean class counts.
#'     }}
#'   \item{\code{variance}}{Variance of the prediction target:
#'     \itemize{
#'       \item For \code{BKP} and \code{TwinBKP} with
#'       \code{type = "probability"}, the posterior variance of the latent
#'       success probability.
#'       \item For \code{BKP} and \code{TwinBKP} with \code{type = "count"},
#'       the posterior predictive variance of a future success count.
#'       \item For \code{DKP} and \code{TwinDKP} with
#'       \code{type = "probability"}, a matrix of marginal posterior variances
#'       for class probabilities.
#'       \item For \code{DKP} and \code{TwinDKP} with \code{type = "count"},
#'       a matrix of marginal posterior predictive variances for class counts.
#'     }}
#'   \item{\code{lower}}{Lower bound of the credible or predictive interval:
#'     \itemize{
#'       \item For \code{BKP} and \code{TwinBKP} with
#'       \code{type = "probability"}, the lower Beta posterior quantile for
#'       the latent success probability.
#'       \item For \code{BKP} and \code{TwinBKP} with \code{type = "count"},
#'       the lower Beta-Binomial posterior predictive quantile for a future
#'       success count.
#'       \item For \code{DKP} and \code{TwinDKP} with
#'       \code{type = "probability"}, a matrix of lower marginal posterior
#'       quantiles for class probabilities.
#'       \item For \code{DKP} and \code{TwinDKP} with \code{type = "count"},
#'       a matrix of lower marginal Beta-Binomial posterior predictive quantiles
#'       for class counts.
#'     }}
#'   \item{\code{upper}}{Upper bound of the credible or predictive interval:
#'     \itemize{
#'       \item For \code{BKP} and \code{TwinBKP} with
#'       \code{type = "probability"}, the upper Beta posterior quantile for
#'       the latent success probability.
#'       \item For \code{BKP} and \code{TwinBKP} with \code{type = "count"},
#'       the upper Beta-Binomial posterior predictive quantile for a future
#'       success count.
#'       \item For \code{DKP} and \code{TwinDKP} with
#'       \code{type = "probability"}, a matrix of upper marginal posterior
#'       quantiles for class probabilities.
#'       \item For \code{DKP} and \code{TwinDKP} with \code{type = "count"},
#'       a matrix of upper marginal Beta-Binomial posterior predictive quantiles
#'       for class counts.
#'     }}
#'   \item{\code{class}}{Predicted class label, returned only when applicable:
#'     \itemize{
#'       \item For \code{BKP} and \code{TwinBKP}, a binary class label based on
#'       the posterior mean and \code{threshold}, returned only for binary
#'       classification data with \code{m = 1} and
#'       \code{type = "probability"}.
#'       \item For \code{DKP} and \code{TwinDKP}, the class with the largest
#'       posterior mean probability, returned only for one-hot classification
#'       data and \code{type = "probability"}.
#'     }}
#'   \item{\code{threshold}}{The classification threshold, returned for
#'     \code{BKP} and \code{TwinBKP} classification predictions.}
#'   \item{\code{CI_level}}{The specified credible interval level.}
#'   \item{\code{type}}{The prediction target, either \code{"probability"} or
#'     \code{"count"}.}
#'   \item{\code{Mnew}}{The trial sizes used for count prediction, returned only
#'     when \code{type = "count"}.}
#'   \item{\code{ess}}{Effective-sample-size calibration method inherited from
#'     the fitted model, returned for \code{BKP} and \code{DKP} objects.}
#'   \item{\code{ess_info}}{ESS diagnostics for \code{BKP} and \code{DKP}
#'     predictions when available. \code{TwinBKP} and \code{TwinDKP} use
#'     uncalibrated global-local posterior updates and do not return ESS
#'     diagnostics.}
#' }
#'
#' @seealso \code{\link{fit_BKP}}, \code{\link{fit_DKP}},
#'   \code{\link{fit_TwinBKP}}, and \code{\link{fit_TwinDKP}} for model fitting;
#'   \code{\link{plot.BKP}}, \code{\link{plot.DKP}},
#'   \code{\link{plot.TwinBKP}}, and \code{\link{plot.TwinDKP}} for visualization.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}. arXiv. <doi:10.48550/arXiv.2508.10447>.
#'
#' @keywords BKP
#'
#' @examples
#' #-------------------------- BKP and TwinBKP ---------------------------
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
#' # A fixed theta is used here only to keep the example fast and reproducible.
#' # In practice, omit theta to select it by leave-one-out cross-validation.
#' model <- fit_BKP(X, y, m, Xbounds = Xbounds, theta = 0.04)
#'
#' # Prediction on training data
#' predict(model)
#'
#' # Prediction on new data
#' Xnew <- matrix(seq(-2, 2, length = 10), ncol = 1) #new data points
#' predict(model, Xnew = Xnew)
#'
#' # Posterior predictive summaries for future success counts
#' Mnew <- sample(100, nrow(Xnew), replace = TRUE)
#' predict(model, Xnew = Xnew, type = "count", Mnew = Mnew)
#'
#' \dontrun{
#' # Larger TwinBKP example
#' n <- 1000
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit TwinBKP model using the default global lengthscale tuning
#' model <- fit_TwinBKP(
#'      X, y, m,
#'      Xbounds = Xbounds
#'    )
#'
#' # Prediction on training data
#' predict(model)
#'
#' # Prediction on new data
#' predict(model, Xnew = Xnew)
#'
#' # Posterior predictive summaries for future success counts
#' predict(model, Xnew = Xnew, type = "count", Mnew = Mnew)
#' }
#'
#' @export
#' @method predict BKP

predict.BKP <- function(object, Xnew = NULL, CI_level = 0.95, threshold = 0.5,
                        type = c("probability", "count"), Mnew = NULL, ...)
{
  # ---- Extract basic information ----
  X <- object$X
  d <- ncol(X)

  # ---- Check and format Xnew ----
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
  n_pred <- if (is.null(Xnew)) nrow(X) else nrow(Xnew)

  # ---- Check scalar probability arguments ----
  if (!is.numeric(CI_level) || length(CI_level) != 1L ||
      is.na(CI_level) || !is.finite(CI_level) ||
      CI_level <= 0 || CI_level >= 1) {
    stop("'CI_level' must be a single finite numeric value strictly between 0 and 1.")
  }

  if (!is.numeric(threshold) || length(threshold) != 1L ||
      is.na(threshold) || !is.finite(threshold) ||
      threshold <= 0 || threshold >= 1) {
    stop("'threshold' must be a single finite numeric value strictly between 0 and 1.")
  }

  type <- match.arg(type)

  # ---- Check Mnew for count prediction ----
  if (type == "count") {
    if (is.null(Mnew)) {
      if (is.null(Xnew)) {
        Mnew <- object$m
      } else {
        stop("When type = 'count' and Xnew is provided, 'Mnew' must also be provided.")
      }
    }

    if (!is.numeric(Mnew)) {
      stop("'Mnew' must be a numeric vector when type = 'count'.")
    }

    if (!(length(Mnew) == 1L || length(Mnew) == n_pred)) {
      stop("'Mnew' must have length 1 or the same length as the number of prediction points.")
    }

    if (anyNA(Mnew) || any(!is.finite(Mnew))) {
      stop("'Mnew' must contain finite values with no NA.")
    }

    if (any(Mnew <= 0) || any(Mnew != floor(Mnew))) {
      stop("'Mnew' must contain positive integers when type = 'count'.")
    }

    Mnew <- as.integer(Mnew)

    if (length(Mnew) == 1L) {
      Mnew <- rep.int(Mnew, n_pred)
    }
  }

  # ---- Posterior parameters ----
  if(!is.null(Xnew)){
    # Extract components
    Xnorm   <- object$Xnorm
    y       <- object$y
    m       <- object$m
    theta   <- object$theta_opt
    kernel  <- object$kernel
    isotropic <- object$isotropic
    prior   <- object$prior
    r0      <- object$r0
    p0      <- object$p0
    Xbounds <- object$Xbounds
    ess     <- if (is.null(object$ess)) "none" else object$ess

    # Normalize Xnew to [0,1]^d
    Xnew_norm <- sweep(Xnew, 2, Xbounds[, 1], "-")
    Xnew_norm <- sweep(Xnew_norm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

    posterior <- bkp_compute_posterior(
      Xquery_norm = Xnew_norm, Xtrain_norm = Xnorm, y = y, m = m,
      theta = theta, kernel = kernel, isotropic = isotropic,
      prior = prior, r0 = r0, p0 = p0, ess = ess
    )
    alpha_n <- posterior$alpha_n
    beta_n <- posterior$beta_n
    ess_info <- posterior$ess_info
  }else{
    # Use stored posterior parameters at training points
    alpha_n <- object$alpha_n
    beta_n  <- object$beta_n
    m <- object$m
    ess_info <- object$ess_info
  }

  # ---- Posterior / posterior predictive summaries ----
  s_n <- alpha_n + beta_n

  if (anyNA(s_n) || any(!is.finite(s_n)) || any(s_n <= 0)) {
    stop("Posterior shape parameters must be positive and finite.")
  }

  # Posterior summaries for the latent success probability pi(x).
  # These quantities are always computed first, because both probability-scale
  # and count-scale predictions are derived from pi(x) | D_n.
  prod_mean <- alpha_n / s_n
  prod_var  <- prod_mean * (1 - prod_mean) / (s_n + 1)

  if (type == "probability") {
    # Prediction target: latent success probability pi(x).
    pred_mean <- prod_mean
    pred_var  <- prod_var

    # Credible intervals for pi(x) | D_n.
    pred_lower <- suppressWarnings(qbeta((1 - CI_level) / 2, alpha_n, beta_n))
    pred_upper <- suppressWarnings(qbeta((1 + CI_level) / 2, alpha_n, beta_n))
  } else {
    # Prediction target: future success count y(x) given Mnew.
    # Marginally, y(x) | D_n, Mnew follows a Beta-Binomial distribution.
    pred_mean <- Mnew * prod_mean
    pred_var <- Mnew * (s_n + Mnew) * prod_var

    # Predictive intervals for y(x) | D_n, Mnew.
    pred_lower <- qbetabinom_rcpp((1 - CI_level) / 2, Mnew, alpha_n, beta_n)
    pred_upper <- qbetabinom_rcpp((1 + CI_level) / 2, Mnew, alpha_n, beta_n)
  }

  # ---- Return object ----
  prediction <- list(
    X = X,
    Xnew = Xnew,
    alpha_n = alpha_n,
    beta_n = beta_n,
    mean     = pred_mean,
    variance = pred_var,
    lower    = pred_lower,
    upper    = pred_upper,
    CI_level = CI_level,
    type = type,
    ess = if (is.null(object$ess)) "none" else object$ess,
    ess_info = ess_info
  )

  if (type == "count") {
    prediction$Mnew <- Mnew
  }

  # Posterior classification label (only for classification data)
  if (type == "probability" && all(m == 1)) {
    prediction$class <- ifelse(pred_mean > threshold, 1, 0)
    prediction$threshold <- threshold
  }

  class(prediction) <- "predict_BKP"
  return(prediction)
}
