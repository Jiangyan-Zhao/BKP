#' @name predict.TwinBKP
#'
#' @title Predict from a Fitted TwinBKP Model
#'
#' @description
#' Prediction method for TwinBKP. For each prediction point, the function runs:
#'
#' 1. **Local point selection**: use a kd-tree to find \code{l_nums} nearest neighbors
#'    from the training set as local points.
#' 2. **Wendland kernel hyperparameter**: compute the coverage radius
#'    \eqn{\hat{\theta}_l} using Equation (14) as the single local-kernel
#'    hyperparameter.
#' 3. **Validation set**: randomly sample \code{v_nums = 2 * g_nums} points from non-global points.
#' 4. **Mixing weight \eqn{\lambda}**: minimize the loss of the mixed kernel
#'    \eqn{\lambda K_g + (1-\lambda) K_l} on the validation set to obtain the
#'    optimal \eqn{\lambda}.
#' 5. **Posterior prediction**: perform BKP posterior updating with the mixed kernel
#'    and return posterior mean/variance/confidence interval.
#'
#' @param object A \code{"TwinBKP"} object returned by \code{fit_TwinBKP()}.
#' @param Xnew Prediction input matrix (unnormalized), with dimension \code{n_new x d}.
#' @param l_nums Positive integer. Number of local neighbors per prediction point.
#'   If \code{NULL} (default), it is set to \code{max(25, 3 * d)}, where \code{d}
#'   is the input dimensionality.
#' @param v_nums Positive integer. Validation set size. If \code{NULL} (default),
#'   it is set to \code{2 * object$g_nums}.
#' @param CI_level Numeric confidence level. Default is \code{0.95}.
#' @param threshold Classification threshold (only used for classification with
#'   \code{m = 1}). Default is \code{0.5}.
#' @param ... Unused.
#'
#' @return A list of class \code{"predict_TwinBKP"} containing:
#' \describe{
#'   \item{\code{X}}{Original training input matrix.}
#'   \item{\code{mean}}{Posterior predictive mean matrix of size \code{n_new x 1}.}
#'   \item{\code{variance}}{Posterior predictive variance matrix of size \code{n_new x 1}.}
#'   \item{\code{lower, upper}}{Lower/upper confidence interval matrices of size \code{n_new x 1}.}
#'   \item{\code{alpha_n, beta_n}}{Posterior Beta shape parameters at prediction points.}
#'   \item{\code{lambda}}{Mixing weight \eqn{\lambda} for each prediction point
#'   (length \code{n_new}).}
#'   \item{\code{theta_l}}{Local kernel hyperparameter for each prediction point
#'   (length \code{n_new}).}
#'   \item{\code{local_idx}}{List of selected local training-row indices for each
#'   prediction point.}
#'   \item{\code{Xnew}}{Original prediction coordinates.}
#'   \item{\code{Xnew_norm}}{Normalized prediction coordinates.}
#' }
#'
#' @seealso \code{\link{fit_TwinBKP}}
#' @examples
#'
#'   # ============================================================== #
#'   # ======================= TwinBKP Examples ====================== #
#'   # ============================================================== #
#'
#'   #-------------------------- 1D Example ---------------------------
#'   set.seed(123)
#'
#'   # Define true success probability function
#'   true_pi_fun <- function(x) {
#'     (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#'   }
#'
#'   n <- 100
#'   Xbounds <- matrix(c(-2, 2), nrow = 1)
#'   X <- tgp::lhs(n = n, rect = Xbounds)
#'   true_pi <- true_pi_fun(X)
#'   m <- sample(100, n, replace = TRUE)
#'   y <- rbinom(n, size = m, prob = true_pi)
#'
#'   # Fit TwinBKP model (global stage only)
#'   model1 <- fit_TwinBKP(X, y, m, Xbounds = Xbounds, g_nums = 10) 
#'
#'   # Prediction on new data
#'   Xnew <- matrix(seq(-2, 2, length.out = 10), ncol = 1)
#'   predict(model1, Xnew = Xnew, l_nums = 10)
#'
#'
#'   #-------------------------- 2D Example ---------------------------
#'   set.seed(123)
#'
#'   # Define 2D latent function and probability transformation
#'   true_pi_fun <- function(X) {
#'     if (is.null(nrow(X))) X <- matrix(X, nrow = 1)
#'     m <- 8.6928
#'     s <- 2.4269
#'     x1 <- 4 * X[, 1] - 2
#'     x2 <- 4 * X[, 2] - 2
#'     a <- 1 + (x1 + x2 + 1)^2 *
#'       (19 - 14 * x1 + 3 * x1^2 - 14 * x2 + 6 * x1 * x2 + 3 * x2^2)
#'     b <- 30 + (2 * x1 - 3 * x2)^2 *
#'       (18 - 32 * x1 + 12 * x1^2 + 48 * x2 - 36 * x1 * x2 + 27 * x2^2)
#'     f <- log(a * b)
#'     f <- (f - m) / s
#'     pnorm(f)
#'   }
#'
#'   n <- 50
#'   Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#'   X <- tgp::lhs(n = n, rect = Xbounds)
#'   true_pi <- true_pi_fun(X)
#'   m <- sample(100, n, replace = TRUE)
#'   y <- rbinom(n, size = m, prob = true_pi)
#'
#'   # Fit TwinBKP model
#'   model2 <- fit_TwinBKP(X, y, m, Xbounds = Xbounds, g_nums = 12) 
#'
#'   # Prediction on new data
#'   x1 <- seq(Xbounds[1, 1], Xbounds[1, 2], length.out = 8)
#'   x2 <- seq(Xbounds[2, 1], Xbounds[2, 2], length.out = 8)
#'   Xnew <- expand.grid(x1 = x1, x2 = x2)
#'   predict(model2, Xnew = Xnew, l_nums = 12)
#' 
#' @export
predict.TwinBKP <- function(
    object,
    Xnew,
    l_nums    = NULL,
    v_nums    = NULL,
    CI_level  = 0.95,
    threshold = 0.5,
    return_type = c("probability", "count"),
    n_trials = NULL,
    ...
) {


  Xnorm        <- object$Xnorm
  y            <- object$y
  m_train      <- object$m
  Xbounds      <- object$Xbounds
  global_idx   <- object$global_idx
  Xnorm_global <- object$Xnorm_global
  y_global     <- object$y_global
  m_global     <- object$m_global
  theta_global <- object$theta_global
  kernel       <- object$kernel
  isotropic    <- object$isotropic
  prior        <- object$prior
  r0           <- object$r0
  p0           <- object$p0
  loss_type    <- object$loss
  g_nums       <- object$g_nums

  n <- nrow(Xnorm)
  d <- ncol(Xnorm)
 

  if (is.null(nrow(Xnew))) Xnew <- matrix(Xnew, nrow = 1)
  Xnew <- as.matrix(Xnew)
  if (!is.numeric(Xnew)) stop("'Xnew' must be numeric.")
  if (ncol(Xnew) != d)   stop("'Xnew' must have the same number of columns as training 'X'.")
  if (anyNA(Xnew))       stop("'Xnew' contains NA values.")

  if (!is.numeric(CI_level) || CI_level <= 0 || CI_level >= 1) {
    stop("'CI_level' must be strictly between 0 and 1.")
  }
  if (!is.numeric(threshold) || threshold <= 0 || threshold >= 1) {
    stop("'threshold' must be strictly between 0 and 1.")
  }
  return_type <- match.arg(return_type)
  if (return_type == "count") {
    if (is.null(n_trials)) {
      stop("When return_type = 'count', 'n_trials' must be provided.")
    }
    if (!is.numeric(n_trials) || length(n_trials) != 1 || n_trials <= 0 || n_trials != as.integer(n_trials)) {
      stop("'n_trials' must be a positive integer when return_type = 'count'.")
    }
    n_trials <- as.integer(n_trials)
  }

  Xnew_norm <- sweep(Xnew, 2, Xbounds[, 1], "-")
  Xnew_norm <- sweep(Xnew_norm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

  n_new <- nrow(Xnew_norm)


  if (is.null(l_nums)) l_nums <- max(25L, 3L * d)
  if (!is.numeric(l_nums) || length(l_nums) != 1 || l_nums <= 0) {
    stop("'l_nums' must be a positive integer.")
  }
  l_nums <- as.integer(l_nums)
  l_eff <- min(l_nums, n)

  if (is.null(g_nums) && !is.null(object$g)) g_nums <- object$g
  if (is.null(g_nums)) {
    stop("Cannot determine 'g_nums' from model object. Please refit with updated 'fit_TwinBKP()' or provide 'v_nums' explicitly.")
  }
  if (is.null(v_nums)) v_nums <- 2L * as.integer(g_nums)
  if (!is.numeric(v_nums) || length(v_nums) != 1 || v_nums <= 0) {
    stop("'v_nums' must be a positive integer.")
  }
  v_nums <- as.integer(v_nums)

  knn_result <- get_knnx_nanoflann_rcpp(
    data = Xnorm,
    query = Xnew_norm,
    k = l_eff
  )
  local_idx_mat  <- knn_result$nn.index
  local_dist_mat <- knn_result$nn.dist

  knn_global <- get_knnx_nanoflann_rcpp(
    data = Xnorm_global,
    query = Xnorm,
    k = 1L
  )
  theta_l <- max(knn_global$nn.dist)

  q_wend <- floor(d / 2) + 3

  non_global_idx <- setdiff(seq_len(n), global_idx)

  v_eff <- min(v_nums, length(non_global_idx))
  if (v_eff < v_nums) {
    warning(sprintf(
      "Only %d non-global points available; validation set size reduced from %d to %d.",
      length(non_global_idx), v_nums, v_eff
    ))
  }

  val_idx   <- sample(non_global_idx, size = v_eff, replace = FALSE)
  Xnorm_val <- Xnorm[val_idx, , drop = FALSE]
  y_val     <- y[val_idx, , drop = FALSE]
  m_val     <- m_train[val_idx, , drop = FALSE]


  K_g_val <- kernel_matrix(
    X = Xnorm_val,
    theta = theta_global,
    kernel = kernel,
    isotropic = isotropic
  )

  pred_mean     <- numeric(n_new)
  pred_var      <- numeric(n_new)
  pred_lower    <- numeric(n_new)
  pred_upper    <- numeric(n_new)
  pred_alpha_n  <- numeric(n_new)
  pred_beta_n   <- numeric(n_new)
  pred_lambda   <- numeric(n_new)
  pred_theta_l  <- rep(theta_l, n_new)
  local_idx_out <- vector("list", n_new)

  for (i in seq_len(n_new)) {

    loc_idx  <- local_idx_mat[i, ]
    Xnorm_loc <- Xnorm[loc_idx, , drop = FALSE]
    y_loc     <- y[loc_idx, , drop = FALSE]
    m_loc     <- m_train[loc_idx, , drop = FALSE]
    local_idx_out[[i]] <- loc_idx

    K_l_val <- wendland_kernel_rcpp(
      X1 = Xnorm_val,
      X2 = Xnorm_val,
      theta = theta_l,
      q_wend = q_wend
    )

    prior_val <- get_prior(
      prior = prior, model = "BKP",
      r0 = r0, p0 = p0,
      y = y_val, m = m_val,
      K = K_g_val
    )
    alpha0_val <- prior_val$alpha0
    beta0_val  <- prior_val$beta0

    opt_lambda <- optimize_lambda_bkp_rcpp(
      K_g = K_g_val,
      K_l = K_l_val,
      y = as.numeric(y_val),
      m = as.numeric(m_val),
      alpha0 = as.numeric(alpha0_val),
      beta0 = as.numeric(beta0_val),
      loss = loss_type
    )
    lambda_i <- as.numeric(opt_lambda$lambda_opt)
    pred_lambda[i] <- lambda_i

    K_g_star <- kernel_matrix(
      X        = matrix(Xnew_norm[i, ], nrow = 1),
      Xprime   = Xnorm_global,
      theta    = theta_global,
      kernel   = kernel,
      isotropic = isotropic
    )  # (1 x g)

    K_l_star <- wendland_kernel_rcpp(
      X1    = matrix(Xnew_norm[i, ], nrow = 1),
      X2    = Xnorm_loc,
      theta = theta_l,
      q_wend = q_wend
    )  # (1 x l)
    prior_g_star <- get_prior(
      prior = prior, model = "BKP",
      r0 = r0, p0 = p0,
      y = y_global, m = m_global,
      K = K_g_star
    )

    prior_l_star <- get_prior(
      prior = prior, model = "BKP",
      r0 = r0, p0 = p0,
      y = y_loc, m = m_loc,
      K = K_l_star
    )

    # alpha_n = alpha0 + lambda * K_g * y_g + (1-lambda) * K_l * y_l
    alpha_n_i <- lambda_i * (as.numeric(prior_g_star$alpha0) +
                               as.numeric(K_g_star %*% as.numeric(y_global))) +
      (1 - lambda_i) * (as.numeric(prior_l_star$alpha0) +
                          as.numeric(K_l_star %*% as.numeric(y_loc)))

    beta_n_i  <- lambda_i * (as.numeric(prior_g_star$beta0) +
                               as.numeric(K_g_star %*% as.numeric(m_global - y_global))) +
      (1 - lambda_i) * (as.numeric(prior_l_star$beta0) +
                          as.numeric(K_l_star %*% as.numeric(m_loc - y_loc)))

    pred_alpha_n[i] <- alpha_n_i
    pred_beta_n[i]  <- beta_n_i

    eps <- 1e-10
    ab_sum <- max(alpha_n_i + beta_n_i, eps)

    if (return_type == "probability") {
      pred_mean[i] <- alpha_n_i / ab_sum
      pred_mean[i] <- min(max(pred_mean[i], eps), 1 - eps)
      pred_var[i]  <- pred_mean[i] * (1 - pred_mean[i]) / (ab_sum + 1)

      pred_lower[i] <- suppressWarnings(qbeta((1 - CI_level) / 2, alpha_n_i, beta_n_i))
      pred_upper[i] <- suppressWarnings(qbeta((1 + CI_level) / 2, alpha_n_i, beta_n_i))
    } else {
      pred_mean[i] <- n_trials * alpha_n_i / ab_sum
      pred_var[i] <- n_trials * alpha_n_i * beta_n_i * (ab_sum + n_trials) /
        (pmax(ab_sum^2, eps) * pmax(ab_sum + 1, eps))
      pred_lower[i] <- betabinom_quantile((1 - CI_level) / 2, n_trials, alpha_n_i, beta_n_i)
      pred_upper[i] <- betabinom_quantile((1 + CI_level) / 2, n_trials, alpha_n_i, beta_n_i)
    }
  }

  mean_mat  <- matrix(as.numeric(pred_mean), ncol = 1)
  var_mat   <- matrix(as.numeric(pred_var), ncol = 1)
  lower_mat <- matrix(as.numeric(pred_lower), ncol = 1)
  upper_mat <- matrix(as.numeric(pred_upper), ncol = 1)

  result <- list(
    X         = object$X,
    Xnew      = Xnew,          
    Xnew_norm = Xnew_norm,   
    alpha_n   = as.numeric(pred_alpha_n),
    beta_n    = as.numeric(pred_beta_n),
    mean      = mean_mat,
    variance  = var_mat,
    lower     = lower_mat,
    upper     = upper_mat,
    lambda    = pred_lambda, 
    theta_l   = pred_theta_l,
    theta_global = theta_global,
    local_idx = local_idx_out, 
    CI_level  = CI_level,
    return_type = return_type,
    l_nums    = l_eff,
    v_nums    = v_eff,
    g_nums    = as.integer(g_nums)
  )

  if (return_type == "count") {
    result$n_trials <- n_trials
  }

  if (return_type == "probability" && all(m_train == 1)) {
    result$class     <- ifelse(as.numeric(mean_mat) > threshold, 1, 0)
    result$threshold <- threshold
  }

  class(result) <- "predict_TwinBKP"
  return(result)
}