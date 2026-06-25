#' @name predict
#'
#' @title Posterior Prediction for BKP, DKP, or TwinBKP Models
#'
#' @description Generates posterior summaries from fitted Beta Kernel Process
#'   (BKP), Dirichlet Kernel Process (DKP), or TwinBKP models at training or
#'   new input locations. For BKP and TwinBKP models, summaries can be returned
#'   either for the latent success probability or for a future success count
#'   under the Beta-Binomial posterior predictive distribution. For DKP models,
#'   summaries can be returned either for the latent class probability vector or
#'   for future class counts using marginal Beta-Binomial posterior predictive
#'   distributions.
#'
#' @param object An object of class \code{"BKP"}, \code{"DKP"}, or
#'   \code{"TwinBKP"}, typically returned by \code{\link{fit_BKP}},
#'   \code{\link{fit_DKP}}, or \code{\link{fit_TwinBKP}}.
#' @param Xnew A numeric vector or matrix of new input locations at which to
#'   generate predictions. If \code{NULL}, predictions are returned for the
#'   training data.
#' @param CI_level Numeric between 0 and 1 specifying the credible level for
#'   posterior intervals (default \code{0.95} for 95% credible interval).
#' @param threshold Numeric between 0 and 1 specifying the classification
#'   threshold for binary predictions based on posterior mean (used only for
#'   BKP; default is \code{0.5}).
#' @param type Character string specifying the prediction target. The default
#'   \code{"probability"} returns posterior summaries for latent success
#'   probabilities (BKP) or latent class probabilities (DKP). The option
#'   \code{"count"} returns posterior predictive summaries for future counts
#'   under the Beta-Binomial distribution (BKP) or marginal Beta-Binomial
#'   distributions for each multinomial class (DKP).
#' @param Mnew Positive integer trial size used when \code{type = "count"}. It
#'   can be either a scalar, applied to all prediction points, or a vector with
#'   the same length as the number of prediction points. If \code{Xnew = NULL}
#'   and \code{Mnew = NULL}, the training trial sizes \code{object$m} (BKP) or
#'   \code{rowSums(object$Y)} (DKP) are used.
#' @param ... Additional arguments passed to generic \code{predict} methods
#'   (currently not used; included for S3 method consistency).
#'
#' @return A list containing posterior or posterior predictive summaries:
#' \describe{
#'   \item{\code{X}}{The original training input locations.}
#'   \item{\code{Xnew}}{The new input locations for prediction. If \code{NULL},
#'     predictions are returned at the training input locations.}
#'   \item{\code{alpha_n}, \code{beta_n}}{Posterior shape parameters:
#'     \itemize{
#'       \item BKP and TwinBKP: Vectors \code{alpha_n} and \code{beta_n} of Beta posterior
#'       shape parameters.
#'       \item DKP: Matrix \code{alpha_n} of Dirichlet posterior concentration
#'       parameters, with rows corresponding to input locations and columns to
#'       classes.
#'     }}
#'   \item{\code{mean}}{Mean of the prediction target:
#'     \itemize{
#'       \item BKP and TwinBKP with \code{type = "probability"}: posterior mean of the latent
#'       success probability.
#'       \item BKP and TwinBKP with \code{type = "count"}: posterior predictive mean of the
#'       future success count.
#'       \item DKP with \code{type = "probability"}: matrix of posterior mean class probabilities.
#'       \item DKP with \code{type = "count"}: matrix of posterior predictive mean class counts.
#'     }}
#'   \item{\code{variance}}{Variance of the prediction target:
#'     \itemize{
#'       \item BKP and TwinBKP with \code{type = "probability"}: posterior variance of the
#'       latent success probability.
#'       \item BKP and TwinBKP with \code{type = "count"}: posterior predictive variance of
#'       the future success count.
#'       \item DKP with \code{type = "probability"}: matrix of posterior variances for class probabilities.
#'       \item DKP with \code{type = "count"}: matrix of posterior predictive variances for class counts.
#'     }}
#'   \item{\code{lower}}{Lower bound of the credible interval:
#'     \itemize{
#'       \item BKP and TwinBKP with \code{type = "probability"}: lower Beta posterior
#'       quantile for the latent success probability.
#'       \item BKP and TwinBKP with \code{type = "count"}: lower Beta-Binomial posterior
#'       predictive quantile for the future success count.
#'       \item DKP with \code{type = "probability"}: matrix of lower posterior quantiles for class probabilities.
#'       \item DKP with \code{type = "count"}: matrix of lower posterior predictive quantiles for class counts.
#'     }}
#'   \item{\code{upper}}{Upper bound of the credible interval:
#'     \itemize{
#'       \item BKP and TwinBKP with \code{type = "probability"}: upper Beta posterior
#'       quantile for the latent success probability.
#'       \item BKP and TwinBKP with \code{type = "count"}: upper Beta-Binomial posterior
#'       predictive quantile for the future success count.
#'       \item DKP with \code{type = "probability"}: matrix of upper posterior quantiles for class probabilities.
#'       \item DKP with \code{type = "count"}: matrix of upper posterior predictive quantiles for class counts.
#'     }}
#'   \item{\code{class}}{Predicted label, where applicable:
#'     \itemize{
#'       \item BKP and TwinBKP: binary class label based on posterior mean and
#'       \code{threshold}, only for binary data with \code{m = 1} and
#'       \code{type = "probability"}.
#'       \item DKP: predicted class label with the highest posterior mean
#'       probability.
#'     }}
#'   \item{\code{CI_level}}{The specified credible interval level.}
#'   \item{\code{type}}{The prediction target used.}
#'   \item{\code{ess}}{Effective-sample-size calibration method inherited from
#'     the fitted model, returned for BKP and DKP objects when available.}
#'   \item{\code{ess_info}}{ESS diagnostics for BKP and DKP predictions when
#'     available. TwinBKP uses the uncalibrated global-local posterior update and
#'     does not return ESS diagnostics.}
#'   \item{\code{Mnew}}{The trial sizes used for count prediction, returned only
#'     when \code{type = "count"}.}
#' }
#'
#' @seealso \code{\link{fit_BKP}}, \code{\link{fit_DKP}}, and
#'   \code{\link{fit_TwinBKP}} for model fitting; \code{\link{plot.BKP}},
#'   \code{\link{plot.DKP}}, and \code{\link{plot.TwinBKP}} for visualization.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}. arXiv. \doi{10.48550/arXiv.2508.10447}
#'
#' @keywords BKP
#'
#' @examples
#' # ============================================================== #
#' # ========================= BKP Examples ======================= #
#' # ============================================================== #
#'
#' #-------------------------- 1D Example ---------------------------
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
#' model1 <- fit_BKP(X, y, m, Xbounds=Xbounds)
#'
#' # Prediction on training data
#' predict(model1)
#'
#' # Prediction on new data
#' Xnew <- matrix(seq(-2, 2, length = 10), ncol=1) #new data points
#' predict(model1, Xnew = Xnew)
#'
#' # Posterior predictive summaries for future success counts
#' Mnew <- sample(100, nrow(Xnew), replace = TRUE)
#' predict(model1, Xnew = Xnew, type = "count", Mnew = Mnew)
#'
#' #-------------------------- 2D Example ---------------------------
#' set.seed(123)
#'
#' # Define 2D latent function and probability transformation
#' true_pi_fun <- function(X) {
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
#'   return(pnorm(f))  # Transform to probability
#' }
#'
#' n <- 100
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit BKP model
#' model2 <- fit_BKP(X, y, m, Xbounds=Xbounds)
#'
#' # Prediction on training data
#' predict(model2)
#'
#' # Prediction on new data
#' x1 <- seq(Xbounds[1,1], Xbounds[1,2], length.out = 10)
#' x2 <- seq(Xbounds[2,1], Xbounds[2,2], length.out = 10)
#' Xnew <- expand.grid(x1 = x1, x2 = x2)
#' predict(model2, Xnew = Xnew)
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

    if (ncol(Xnew) != d) {
      stop("The number of columns in 'Xnew' must match the original input dimension.")
    }
  }
  n_pred <- if (is.null(Xnew)) nrow(X) else nrow(Xnew)

  # ---- Check scalar probability arguments ----
  if (!is.numeric(CI_level) || length(CI_level) != 1 || CI_level <= 0 || CI_level >= 1) {
    stop("'CI_level' must be a single numeric value strictly between 0 and 1.")
  }

  if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0 || threshold >= 1) {
    stop("'threshold' must be a single numeric value strictly between 0 and 1.")
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

    posterior <- .bkp_compute_posterior(
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
