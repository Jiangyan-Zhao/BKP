#' @rdname predict
#' @keywords TwinDKP
#' @export
#' @method predict TwinDKP
predict.TwinDKP <- function(object, Xnew = NULL, CI_level = 0.95,
                            type = c("probability", "count"),
                            Mnew = NULL, ...) {
  X <- object$X
  d <- ncol(X)
  type <- match.arg(type)

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
    if (anyNA(Xnew) || any(!is.finite(Xnew))) {
      stop("'Xnew' must contain only finite values with no NA, NaN, or Inf.")
    }
  }

  n_pred <- if (is.null(Xnew)) nrow(X) else nrow(Xnew)

  if (!is.numeric(CI_level) || length(CI_level) != 1 ||
      CI_level <= 0 || CI_level >= 1) {
    stop("'CI_level' must be a single numeric value strictly between 0 and 1.")
  }

  # ---- Check Mnew for count prediction ----
  if (type == "count") {
    if (is.null(Mnew)) {
      if (is.null(Xnew)) {
        Mnew <- rowSums(object$Y)
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

    Mnew <- as.integer(if (length(Mnew) == 1L) rep.int(Mnew, n_pred) else Mnew)
  }

  # ---- Posterior parameters ----
  if (is.null(Xnew)) {
    alpha_n <- object$alpha_n
  } else {
    Xnew_norm <- sweep(Xnew, 2, object$Xbounds[, 1], "-")
    Xnew_norm <- sweep(Xnew_norm, 2, object$Xbounds[, 2] - object$Xbounds[, 1], "/")

    leaf_size <- if (is.null(object$control$leaf_size)) 8L else object$control$leaf_size
    local_indices <- twin_local_indices_rcpp(
      Xtrain_norm = object$Xnorm,
      Xquery_norm = Xnew_norm,
      g_indices = object$global_indices,
      l = object$control$l,
      leaf_size = leaf_size
    )

    posterior <- twin_dkp_posterior_rcpp(
      Xquery_norm = Xnew_norm,
      Xtrain_norm = object$Xnorm,
      Y = object$Y,
      g_indices = as.integer(object$global_indices),
      local_indices = local_indices,
      theta_g = as.numeric(object$theta_g),
      theta_l = as.numeric(object$theta_l),
      global_kernel = object$global_kernel,
      local_kernel = object$local_kernel,
      isotropic = isTRUE(object$isotropic),
      prior = object$prior,
      r0 = object$r0,
      p0 = as.numeric(object$p0),
      store_kernel = FALSE
    )
    alpha_n <- posterior$alpha_n
  }

  # ---- Posterior / posterior predictive summaries ----
  alpha_n <- as.matrix(alpha_n)
  if (anyNA(alpha_n) || any(!is.finite(alpha_n)) || any(alpha_n <= 0)) {
    stop("Posterior concentration parameters 'alpha_n' must be positive and finite.")
  }

  A <- rowSums(alpha_n)
  if (anyNA(A) || any(!is.finite(A)) || any(A <= 0)) {
    stop("Posterior concentration row sums must be positive and finite.")
  }

  Amat <- matrix(A, nrow(alpha_n), ncol(alpha_n))
  beta_n <- Amat - alpha_n
  prob_mean <- alpha_n / Amat
  prob_var <- prob_mean * (1 - prob_mean) / (Amat + 1)

  class_names <- paste0("class", seq_len(ncol(alpha_n)))

  if (type == "probability") {
    pred_mean <- prob_mean
    pred_var <- prob_var
    pred_lower <- qbeta((1 - CI_level) / 2, alpha_n, beta_n)
    pred_upper <- qbeta((1 + CI_level) / 2, alpha_n, beta_n)
  } else {
    Mmat <- matrix(Mnew, nrow(alpha_n), ncol(alpha_n))
    pred_mean <- Mmat * prob_mean
    pred_var <- Mmat * (Amat + Mmat) * prob_var
    pred_lower <- matrix(
      qbetabinom_rcpp(
        (1 - CI_level) / 2,
        as.vector(Mmat),
        as.vector(alpha_n),
        as.vector(beta_n)
      ),
      nrow(alpha_n),
      ncol(alpha_n)
    )
    pred_upper <- matrix(
      qbetabinom_rcpp(
        (1 + CI_level) / 2,
        as.vector(Mmat),
        as.vector(alpha_n),
        as.vector(beta_n)
      ),
      nrow(alpha_n),
      ncol(alpha_n)
    )
  }

  colnames(pred_mean) <- class_names
  colnames(pred_var) <- class_names
  colnames(pred_lower) <- class_names
  colnames(pred_upper) <- class_names

  out <- list(
    X = object$X,
    Xnew = Xnew,
    alpha_n = alpha_n,
    mean = pred_mean,
    variance = pred_var,
    lower = pred_lower,
    upper = pred_upper,
    CI_level = CI_level,
    type = type,
    ess = "none",
    ess_info = NULL
  )

  if (type == "count") {
    out$Mnew <- Mnew
  }
  class(out) <- "predict_TwinDKP"
  out
}
