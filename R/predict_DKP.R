#' @rdname predict
#'
#' @keywords DKP
#'
#' @examples
#' #-------------------------- DKP and TwinDKP ---------------------------
#' # Define true class probability function (3-class)
#' true_pi_fun <- function(X) {
#'   p1 <- 1/(1+exp(-3*X))
#'   p2 <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
#'   return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
#' }
#'
#' n <- 30
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(150, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit DKP model
#' model <- fit_DKP(X, Y, Xbounds = Xbounds)
#'
#' # Prediction on training data
#' predict(model)
#'
#' # Prediction on new data
#' Xnew = matrix(seq(-2, 2, length = 10), ncol=1) #new data points
#' predict(model, Xnew)
#'
#' \dontrun{
#' # Larger TwinDKP example
#' n <- 1000
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(150, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit TwinDKP model
#' model <- fit_TwinDKP(X, Y, Xbounds = Xbounds)
#'
#' # Prediction on training data
#' predict(model)
#'
#' # Prediction on new data
#' predict(model, Xnew)
#' }
#'
#' @export
#' @method predict DKP

predict.DKP <- function(object, Xnew = NULL, CI_level = 0.95,
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

  type <- match.arg(type)

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

    Mnew <- as.integer(Mnew)

    if (length(Mnew) == 1L) {
      Mnew <- rep.int(Mnew, n_pred)
    }
  }

  # ---- Posterior parameters ----
  if(!is.null(Xnew)){
    # Extract components
    Xnorm   <- object$Xnorm
    Y       <- object$Y
    theta   <- object$theta_opt
    ess_info <- object$ess_info
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

    posterior <- .dkp_compute_posterior(
      Xquery_norm = Xnew_norm, Xtrain_norm = Xnorm, Y = Y, theta = theta,
      kernel = kernel, isotropic = isotropic, prior = prior, r0 = r0,
      p0 = p0, ess = ess
    )
    alpha_n <- posterior$alpha_n
    ess_info <- posterior$ess_info
  }else{
    # Use stored posterior parameters at training points
    alpha_n <- object$alpha_n
    Y       <- object$Y
    ess_info <- object$ess_info
  }

  # ---- Posterior / posterior predictive summaries ----
  alpha_n <- as.matrix(alpha_n)

  if (anyNA(alpha_n) || any(!is.finite(alpha_n)) || any(alpha_n <= 0)) {
    stop("Posterior concentration parameters 'alpha_n' must be positive and finite.")
  }

  # Total Dirichlet concentration parameter A_n(x) for each prediction point.
  row_sum <- rowSums(alpha_n)

  # Matrix form of A_n(x), used for element-wise operations over classes.
  A_mat <- matrix(row_sum, nrow = nrow(alpha_n), ncol = ncol(alpha_n))

  # For each class s, the marginal posterior distribution of pi_s(x) is
  # Beta(alpha_{n,s}(x), A_n(x) - alpha_{n,s}(x)).
  beta_n <- A_mat - alpha_n

  # Posterior summaries for latent class probabilities pi_s(x).
  # These are probability-scale quantities and are used as building blocks for
  # both probability-scale and count-scale predictions.
  prob_mean <- alpha_n / A_mat
  prob_var  <- prob_mean * (1 - prob_mean) / (A_mat + 1)

  if (type == "probability") {
    # Prediction target: latent class probability pi_s(x).
    # The returned mean, variance, and credible intervals are on the probability scale.
    pred_mean <- prob_mean
    pred_var  <- prob_var

    # Marginal posterior credible intervals for each class probability pi_s(x).
    pred_lower <- matrix(
      suppressWarnings(qbeta((1 - CI_level) / 2, alpha_n, beta_n)),
      nrow = nrow(alpha_n),
      ncol = ncol(alpha_n)
    )

    pred_upper <- matrix(
      suppressWarnings(qbeta((1 + CI_level) / 2, alpha_n, beta_n)),
      nrow = nrow(alpha_n),
      ncol = ncol(alpha_n)
    )
  } else {
    # Prediction target: future class count y_s(x) given Mnew trials.
    # Marginally, each class count follows a Beta-Binomial distribution:
    # y_s(x) | D_n, Mnew ~ Beta-Binomial(Mnew, alpha_{n,s}, A_n - alpha_{n,s}).
    M_mat <- matrix(Mnew, nrow = nrow(alpha_n), ncol = ncol(alpha_n))

    # Posterior predictive mean and variance of the future class counts.
    pred_mean <- M_mat * prob_mean
    pred_var  <- M_mat * (A_mat + M_mat) * prob_var

    # Marginal posterior predictive intervals for each future class count.
    pred_lower <- matrix(
      qbetabinom_rcpp(
        prob = (1 - CI_level) / 2,
        size = as.vector(M_mat),
        alpha = as.vector(alpha_n),
        beta = as.vector(beta_n)
      ),
      nrow = nrow(alpha_n),
      ncol = ncol(alpha_n)
    )

    pred_upper <- matrix(
      qbetabinom_rcpp(
        prob = (1 + CI_level) / 2,
        size = as.vector(M_mat),
        alpha = as.vector(alpha_n),
        beta = as.vector(beta_n)
      ),
      nrow = nrow(alpha_n),
      ncol = ncol(alpha_n)
    )
  }

  class_names <- paste0("class", seq_len(ncol(alpha_n)))
  colnames(pred_mean) <- class_names
  colnames(pred_var) <- class_names
  colnames(pred_lower) <- class_names
  colnames(pred_upper) <- class_names

  # Return structured output
  prediction <- list(
    X        = X,
    Xnew     = Xnew,
    alpha_n  = alpha_n,
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
  if (type == "probability" && all(rowSums(Y) == 1)) {
    prediction$class <- max.col(pred_mean)
  }

  class(prediction) <- "predict_DKP"
  return(prediction)
}
