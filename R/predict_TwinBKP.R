#' @rdname predict
#'
#' @keywords TwinBKP
#'
#' @examples
#' # ============================================================== #
#' # ======================= TwinBKP Examples ===================== #
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
