#' @title Posterior Prediction for TwinBKP Models
#'
#' @description
#' Generates posterior predictive summaries from a fitted Twin Beta Kernel
#' Process (TwinBKP) model at training or new input locations. For the latent
#' success probability, the method returns Beta posterior summaries, including
#' posterior means, variances, and credible intervals. For future binomial
#' success counts, it returns Beta-Binomial posterior predictive summaries.
#'
#' @param object A fitted object of class \code{"TwinBKP"}, typically returned
#'   by \code{\link{fit_TwinBKP}}.
#' @param Xnew A numeric vector or matrix of new input locations at which to
#'   generate predictions. If \code{NULL}, predictions are returned at the
#'   original training locations using the posterior parameters stored in
#'   \code{object}. For one-dimensional input, a numeric vector is accepted and
#'   internally converted to a one-column matrix.
#' @param CI_level A single numeric value strictly between 0 and 1 specifying
#'   the credible level for posterior or posterior predictive intervals. The
#'   default is \code{0.95}.
#' @param threshold A single numeric value strictly between 0 and 1 specifying
#'   the classification threshold applied to the posterior mean. It is used only
#'   when \code{type = "probability"} and all training trial sizes are one.
#'   The default is \code{0.5}.
#' @param type Character string specifying the prediction target. The default
#'   \code{"probability"} returns posterior summaries for the latent success
#'   probability. The option \code{"count"} returns posterior predictive
#'   summaries for future success counts under the Beta-Binomial distribution.
#' @param Mnew Positive integer trial size used when \code{type = "count"}. It
#'   can be either a scalar, applied to all prediction locations, or a vector
#'   with length equal to the number of prediction locations. If
#'   \code{Xnew = NULL} and \code{Mnew = NULL}, the training trial sizes
#'   \code{object$m} are used. If \code{Xnew} is provided and
#'   \code{type = "count"}, \code{Mnew} must be supplied.
#' @param ... Additional arguments passed to the generic \code{\link{predict}}
#'   method. Currently unused.
#'
#' @details
#' TwinBKP is a scalable global-local approximation to the full BKP model.
#' When predictions are requested at new input locations, the new inputs are
#' first normalized using the bounds stored in the fitted object. The prediction
#' at each new location is then computed using the selected global subset and
#' a set of nearest non-global local neighbours. This preserves the BKP
#' conjugate Beta update while avoiding the full \eqn{n}-point kernel
#' aggregation used by \code{\link{predict.BKP}}.
#'
#' For \code{type = "probability"}, the returned mean and variance correspond
#' to the posterior Beta distribution of the latent success probability
#' \eqn{\pi(x)}. The lower and upper bounds are Beta posterior quantiles.
#'
#' For \code{type = "count"}, the returned mean, variance, and interval bounds
#' correspond to the posterior predictive distribution of a future success
#' count with trial size \code{Mnew}, using the Beta-Binomial distribution.
#'
#' @return
#' An object of class \code{"predict_TwinBKP"}, represented as a list with the
#' following components:
#' \describe{
#'   \item{\code{X}}{The original training input matrix.}
#'   \item{\code{Xnew}}{The prediction input matrix. If \code{NULL},
#'     predictions are evaluated at the training locations.}
#'   \item{\code{alpha_n}, \code{beta_n}}{Posterior Beta shape parameters at
#'     the prediction locations.}
#'   \item{\code{mean}}{Posterior mean of the latent success probability when
#'     \code{type = "probability"}, or posterior predictive mean of the future
#'     success count when \code{type = "count"}.}
#'   \item{\code{variance}}{Posterior variance of the latent success probability
#'     when \code{type = "probability"}, or posterior predictive variance of
#'     the future success count when \code{type = "count"}.}
#'   \item{\code{lower}, \code{upper}}{Lower and upper bounds of the credible
#'     interval or posterior predictive interval.}
#'   \item{\code{CI_level}}{The credible interval level used.}
#'   \item{\code{type}}{The prediction target, either \code{"probability"} or
#'     \code{"count"}.}
#'   \item{\code{Mnew}}{The trial sizes used for count prediction. Returned only
#'     when \code{type = "count"}.}
#'   \item{\code{class}}{Predicted binary class labels. Returned only when
#'     \code{type = "probability"} and all training trial sizes are one.}
#'   \item{\code{threshold}}{The classification threshold. Returned only when
#'     binary class labels are produced.}
#' }
#'
#' @seealso
#' \code{\link{fit_TwinBKP}} for model fitting,
#' \code{\link{predict.BKP}} for prediction from the full BKP model,
#' \code{\link{plot.TwinBKP}} for visualization, and
#' \code{\link{simulate.BKP}} for posterior simulation from full BKP objects.
#'
#' @keywords TwinBKP
#'
#' @examples
#' #-------------------------- 1D Example ---------------------------
#' set.seed(123)
#'
#' # Define true success probability function
#' true_pi_fun <- function(x) {
#'   (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' }
#'
#' n <- 1000
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit TwinBKP model
#' model1 <- fit_TwinBKP(X, y, m, Xbounds = Xbounds)
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
#' n <- 1000
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit TwinBKP model
#' model2 <- fit_TwinBKP(X, y, m, Xbounds=Xbounds)
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
#' @method predict TwinBKP
predict.TwinBKP <- function(object, Xnew = NULL, CI_level = 0.95,
                            threshold = 0.5,
                            type = c("probability", "count"),
                            Mnew = NULL, ...) {
  X <- object$X
  d <- ncol(X)

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

  if (!is.numeric(CI_level) || length(CI_level) != 1 || CI_level <= 0 || CI_level >= 1) {
    stop("'CI_level' must be a single numeric value strictly between 0 and 1.")
  }

  if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0 || threshold >= 1) {
    stop("'threshold' must be a single numeric value strictly between 0 and 1.")
  }

  type <- match.arg(type)
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

  if (is.null(Xnew)) {
    alpha_n <- object$alpha_n
    beta_n <- object$beta_n
    m <- object$m
  } else {
    Xnew_norm <- sweep(Xnew, 2, object$Xbounds[, 1], "-")
    Xnew_norm <- sweep(Xnew_norm, 2, object$Xbounds[, 2] - object$Xbounds[, 1], "/")

    leaf_size <- if (!is.null(object$control$leaf_size)) object$control$leaf_size else 8L

    local_indices <- twin_local_indices_rcpp(
      Xtrain_norm = object$Xnorm,
      Xquery_norm = Xnew_norm,
      g_indices = object$global_indices,
      l = object$control$l,
      leaf_size = leaf_size
    )

    posterior <- twin_bkp_posterior_rcpp(
      Xquery_norm = Xnew_norm,
      Xtrain_norm = object$Xnorm,
      y = as.numeric(object$y),
      m = as.numeric(object$m),
      g_indices = as.integer(object$global_indices),
      local_indices = local_indices,
      theta_g = as.numeric(object$theta_g),
      theta_l = as.numeric(object$theta_l),
      global_kernel = object$global_kernel,
      local_kernel = object$local_kernel,
      isotropic = isTRUE(object$isotropic),
      prior = object$prior,
      r0 = object$r0,
      p0 = object$p0,
      store_kernel = FALSE
    )

    alpha_n <- posterior$alpha_n
    beta_n <- posterior$beta_n
    m <- object$m
  }

  s_n <- alpha_n + beta_n
  if (anyNA(s_n) || any(!is.finite(s_n)) || any(s_n <= 0)) {
    stop("Posterior shape parameters must be positive and finite.")
  }

  prod_mean <- alpha_n / s_n
  prod_var <- prod_mean * (1 - prod_mean) / (s_n + 1)

  if (type == "probability") {
    pred_mean <- prod_mean
    pred_var <- prod_var
    pred_lower <- suppressWarnings(qbeta((1 - CI_level) / 2, alpha_n, beta_n))
    pred_upper <- suppressWarnings(qbeta((1 + CI_level) / 2, alpha_n, beta_n))
  } else {
    pred_mean <- Mnew * prod_mean
    pred_var <- Mnew * (s_n + Mnew) * prod_var
    pred_lower <- qbetabinom_rcpp((1 - CI_level) / 2, Mnew, alpha_n, beta_n)
    pred_upper <- qbetabinom_rcpp((1 + CI_level) / 2, Mnew, alpha_n, beta_n)
  }

  prediction <- list(
    X = object$X,
    Xnew = Xnew,
    alpha_n = alpha_n,
    beta_n = beta_n,
    mean = pred_mean,
    variance = pred_var,
    lower = pred_lower,
    upper = pred_upper,
    CI_level = CI_level,
    type = type
  )

  if (type == "count") {
    prediction$Mnew <- Mnew
  }

  if (type == "probability" && all(m == 1)) {
    prediction$class <- ifelse(pred_mean > threshold, 1, 0)
    prediction$threshold <- threshold
  }

  class(prediction) <- "predict_TwinBKP"
  prediction
}
