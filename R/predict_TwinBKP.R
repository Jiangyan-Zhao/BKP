#' Predict from a TwinBKP model
#'
#' @param object A fitted object of class \code{"TwinBKP"}.
#' @param Xnew Optional numeric matrix of new input locations.
#' @param CI_level Credible interval level.
#' @param threshold Classification threshold used when all training trial counts are one.
#' @param type Prediction scale: \code{"probability"} or \code{"count"}.
#' @param Mnew Trial counts for count prediction.
#' @param ... Additional arguments, currently unused.
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
    ess_info <- object$ess_info
    m <- object$m
  } else {
    Xnew_norm <- sweep(Xnew, 2, object$Xbounds[, 1], "-")
    Xnew_norm <- sweep(Xnew_norm, 2, object$Xbounds[, 2] - object$Xbounds[, 1], "/")

    local_indices <- twin_local_indices_rcpp(
      Xtrain_norm = object$Xnorm,
      Xquery_norm = Xnew_norm,
      g_indices = object$global_indices,
      l = object$l,
      leaf_size = object$leaf_size
    )

    posterior <- .twin_bkp_compute_posterior(
      Xquery_norm = Xnew_norm,
      Xtrain_norm = object$Xnorm,
      y = object$y,
      m = object$m,
      g_indices = object$global_indices,
      local_indices = local_indices,
      theta_g = object$theta_g,
      theta_l = object$theta_l,
      global_kernel = object$global_kernel,
      local_kernel = object$local_kernel,
      isotropic = object$isotropic,
      prior = object$prior,
      r0 = object$r0,
      p0 = object$p0,
      ess = object$ess,
      store_kernel = FALSE
    )

    alpha_n <- posterior$alpha_n
    beta_n <- posterior$beta_n
    ess_info <- posterior$ess_info
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
    type = type,
    ess = object$ess,
    ess_info = ess_info
  )

  if (type == "count") {
    prediction$Mnew <- Mnew
  }

  if (type == "probability" && all(m == 1)) {
    prediction$class <- ifelse(pred_mean > threshold, 1, 0)
    prediction$threshold <- threshold
  }

  class(prediction) <- c("predict_TwinBKP", "predict_BKP")
  prediction
}
