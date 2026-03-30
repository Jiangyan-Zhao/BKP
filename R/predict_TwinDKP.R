#' @name predict.TwinDKP
#'
#' @title Predict from a Fitted TwinDKP Model
#'
#' @description
#' Prediction method for TwinDKP. For each prediction point, the function runs:
#' \enumerate{
#'   \item Local point selection via kNN with \code{l_nums} neighbors.
#'   \item Local Wendland kernel construction with hyperparameter \code{theta_l}.
#'   \item Validation subset sampling with size \code{v_nums}.
#'   \item Mixing-weight optimization for \eqn{\lambda K_g + (1-\lambda)K_l}.
#'   \item Posterior prediction for class probabilities and uncertainty.
#' }
#'
#' @param object A \code{"TwinDKP"} object returned by \code{fit_TwinDKP()}.
#' @param Xnew Prediction input matrix (unnormalized), with dimension
#'   \code{n_new x d}. If \code{NULL}, use training \code{X}.
#' @param l_nums Positive integer. Number of local neighbors per prediction
#'   point. If \code{NULL} (default), set to \code{max(25, 3 * d)}.
#' @param v_nums Positive integer. Validation set size. If \code{NULL}
#'   (default), set to \code{2 * object$g_nums}.
#' @param CI_level Numeric confidence level in \code{(0,1)}. Default is
#'   \code{0.95}.
#' @param ... Unused.
#'
#' @return A list of class \code{"predict_TwinDKP"} containing:
#' \describe{
#'   \item{\code{X}}{Original training input matrix.}
#'   \item{\code{Xnew}}{Prediction input matrix used.}
#'   \item{\code{Xnew_norm}}{Normalized prediction matrix.}
#'   \item{\code{alpha_n}}{Posterior Dirichlet parameters at \code{Xnew}
#'     (\code{n_new x q}).}
#'   \item{\code{mean}}{Posterior mean class probabilities (\code{n_new x q}).}
#'   \item{\code{variance}}{Posterior variances (\code{n_new x q}).}
#'   \item{\code{lower, upper}}{Marginal credible interval bounds
#'     (\code{n_new x q}).}
#'   \item{\code{lambda}}{Optimized mixing weights (length \code{n_new}).}
#'   \item{\code{theta_l}}{Local kernel hyperparameter (length \code{n_new}).}
#'   \item{\code{local_idx}}{List of local training indices for each point.}
#'   \item{\code{CI_level}}{Credible interval level.}
#'   \item{\code{l_nums, v_nums, g_nums}}{Effective tuning sizes used.}
#'   \item{\code{class}}{Predicted class labels when \code{rowSums(Y)=1}.}
#' }
#'
#' @seealso \code{\link{fit_TwinDKP}}, \code{\link{predict.DKP}}
#' @importFrom stats optimise qbeta
#' @export
predict.TwinDKP <- function(
    object,
    Xnew = NULL,
    l_nums = NULL,
    v_nums = NULL,
    CI_level = 0.95,
    ...
) {
  Xnorm        <- object$Xnorm
  Y_train      <- object$Y
  Xbounds      <- object$Xbounds
  global_idx   <- object$global_idx
  Xnorm_global <- object$Xnorm_global
  Y_global     <- object$Y_global
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
  q <- ncol(Y_train)

  if (is.null(Xnew)) Xnew <- object$X
  if (is.null(nrow(Xnew))) Xnew <- matrix(Xnew, nrow = 1)
  Xnew <- as.matrix(Xnew)
  if (!is.numeric(Xnew)) stop("'Xnew' must be numeric.")
  if (ncol(Xnew) != d) stop("'Xnew' must have the same number of columns as training 'X'.")
  if (anyNA(Xnew)) stop("'Xnew' contains NA values.")

  if (!is.numeric(CI_level) || length(CI_level) != 1 || CI_level <= 0 || CI_level >= 1) {
    stop("'CI_level' must be strictly between 0 and 1.")
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
  local_idx_mat <- knn_result$nn.index

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

  val_idx <- sample(non_global_idx, size = v_eff, replace = FALSE)
  Xnorm_val <- Xnorm[val_idx, , drop = FALSE]
  Y_val <- Y_train[val_idx, , drop = FALSE]

  K_g_val <- kernel_matrix(
    X = Xnorm_val,
    theta = theta_global,
    kernel = kernel,
    isotropic = isotropic
  )

  .wendland_kernel <- function(X1, X2, theta) {
    n1 <- nrow(X1); n2 <- nrow(X2)
    K <- matrix(0, n1, n2)
    for (i in seq_len(n1)) {
      r <- sqrt(rowSums((X2 - matrix(X1[i, ], nrow = n2, ncol = ncol(X1), byrow = TRUE))^2))
      u <- r / theta
      K[i, ] <- (q_wend * u + 1) * pmax(0, 1 - u)^q_wend
    }
    K
  }

  K_l_val <- .wendland_kernel(Xnorm_val, Xnorm_val, theta_l)
  alpha0_val <- get_prior(
    prior = prior, model = "DKP",
    r0 = r0, p0 = p0,
    Y = Y_val, K = K_g_val
  )

  .mixed_loss <- function(lambda, K_g, K_l, Y_v, alpha0_v) {
    lambda <- pmin(pmax(lambda, 0), 1)
    K_mix <- lambda * K_g + (1 - lambda) * K_l
    if (loss_type == "brier") {
      loss_fun_brier_dkp_rcpp(K_mix, as.matrix(Y_v), as.matrix(alpha0_v))
    } else {
      loss_fun_logloss_dkp_rcpp(K_mix, as.matrix(Y_v), as.matrix(alpha0_v))
    }
  }

  pred_alpha_n <- matrix(0, n_new, q)
  pred_mean <- matrix(0, n_new, q)
  pred_var <- matrix(0, n_new, q)
  pred_lower <- matrix(0, n_new, q)
  pred_upper <- matrix(0, n_new, q)
  pred_lambda <- numeric(n_new)
  pred_theta_l <- rep(theta_l, n_new)
  local_idx_out <- vector("list", n_new)

  for (i in seq_len(n_new)) {
    loc_idx <- local_idx_mat[i, ]
    Xnorm_loc <- Xnorm[loc_idx, , drop = FALSE]
    Y_loc <- Y_train[loc_idx, , drop = FALSE]
    local_idx_out[[i]] <- loc_idx

    opt_lambda <- optimise(
      f = .mixed_loss,
      interval = c(0, 1),
      K_g = K_g_val,
      K_l = K_l_val,
      Y_v = Y_val,
      alpha0_v = alpha0_val,
      maximum = FALSE
    )
    lambda_i <- opt_lambda$minimum
    pred_lambda[i] <- lambda_i

    K_g_star <- kernel_matrix(
      X = matrix(Xnew_norm[i, ], nrow = 1),
      Xprime = Xnorm_global,
      theta = theta_global,
      kernel = kernel,
      isotropic = isotropic
    )
    K_l_star <- .wendland_kernel(
      X1 = matrix(Xnew_norm[i, ], nrow = 1),
      X2 = Xnorm_loc,
      theta = theta_l
    )

    alpha0_g_star <- get_prior(
      prior = prior, model = "DKP",
      r0 = r0, p0 = p0,
      Y = Y_global, K = K_g_star
    )

    alpha0_l_star <- get_prior(
      prior = prior, model = "DKP",
      r0 = r0, p0 = p0,
      Y = Y_loc, K = K_l_star
    )

    alpha_n_i <- lambda_i * (as.numeric(alpha0_g_star) + as.numeric(K_g_star %*% Y_global)) +
      (1 - lambda_i) * (as.numeric(alpha0_l_star) + as.numeric(K_l_star %*% Y_loc))

    alpha_n_i <- pmax(alpha_n_i, 1e-10)
    row_sum_i <- sum(alpha_n_i)
    mean_i <- alpha_n_i / pmax(row_sum_i, 1e-10)
    mean_i <- pmin(pmax(mean_i, 1e-10), 1 - 1e-10)
    var_i <- mean_i * (1 - mean_i) / (row_sum_i + 1)
    beta_i <- pmax(row_sum_i - alpha_n_i, 1e-10)

    pred_alpha_n[i, ] <- alpha_n_i
    pred_mean[i, ] <- mean_i
    pred_var[i, ] <- var_i
    pred_lower[i, ] <- suppressWarnings(qbeta((1 - CI_level) / 2, alpha_n_i, beta_i))
    pred_upper[i, ] <- suppressWarnings(qbeta((1 + CI_level) / 2, alpha_n_i, beta_i))
  }

  class_names <- paste0("class", seq_len(q))
  colnames(pred_alpha_n) <- class_names
  colnames(pred_mean) <- class_names
  colnames(pred_var) <- class_names
  colnames(pred_lower) <- class_names
  colnames(pred_upper) <- class_names

  result <- list(
    X = object$X,
    Xnew = Xnew,
    Xnew_norm = Xnew_norm,
    alpha_n = pred_alpha_n,
    mean = pred_mean,
    variance = pred_var,
    lower = pred_lower,
    upper = pred_upper,
    lambda = pred_lambda,
    theta_l = pred_theta_l,
    theta_global = theta_global,
    local_idx = local_idx_out,
    CI_level = CI_level,
    l_nums = l_eff,
    v_nums = v_eff,
    g_nums = as.integer(g_nums)
  )

  if (all(rowSums(Y_train) == 1)) {
    result$class <- max.col(pred_mean)
  }

  class(result) <- "predict_TwinDKP"
  result
}
