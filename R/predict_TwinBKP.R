#' @name predict.TwinBKP
#'
#' @title Predict from a Fitted TwinBKP Model
#'
#' @description
#' Prediction method for TwinBKP. For each prediction point, the function runs:
#'
#' 1. **Local point selection**: use a kd-tree to find \code{l} nearest neighbors
#'    from the training set as local points.
#' 2. **Wendland kernel hyperparameter**: compute the coverage radius
#'    \eqn{\hat{\theta}_l} using Equation (14) as the single local-kernel
#'    hyperparameter.
#' 3. **Validation set**: randomly sample \code{v = 2g} points from non-global points.
#' 4. **Mixing weight \eqn{\lambda}**: minimize the loss of the mixed kernel
#'    \eqn{\lambda K_g + (1-\lambda) K_l} on the validation set to obtain the
#'    optimal \eqn{\lambda}.
#' 5. **Posterior prediction**: perform BKP posterior updating with the mixed kernel
#'    and return posterior mean/variance/confidence interval.
#'
#' @param object A \code{"TwinBKP"} object returned by \code{fit_TwinBKP()}.
#' @param Xnew Prediction input matrix (unnormalized), with dimension \code{n_new x d}.
#' @param l Positive integer. Number of local neighbors per prediction point.
#'   Default is \code{50}.
#' @param v Positive integer. Validation set size. Default is \code{2 * object$g}.
#' @param CI_level Numeric confidence level. Default is \code{0.95}.
#' @param threshold Classification threshold (only used for classification with
#'   \code{m = 1}). Default is \code{0.5}.
#' @param ... Unused.
#'
#' @return A list of class \code{"predict_TwinBKP"} containing:
#' \describe{
#'   \item{\code{mean}}{Posterior predictive mean vector (length \code{n_new}).}
#'   \item{\code{variance}}{Posterior predictive variance vector.}
#'   \item{\code{lower, upper}}{Lower/upper confidence interval vectors.}
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
#' if (requireNamespace("FNN", quietly = TRUE)) {
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
#'   model1 <- fit_TwinBKP(X, y, m, Xbounds = Xbounds, g = 10) 
#'
#'   # Prediction on new data
#'   Xnew <- matrix(seq(-2, 2, length.out = 10), ncol = 1)
#'   predict(model1, Xnew = Xnew, l = 10)
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
#'   model2 <- fit_TwinBKP(X, y, m, Xbounds = Xbounds, g = 12) 
#'
#'   # Prediction on new data
#'   x1 <- seq(Xbounds[1, 1], Xbounds[1, 2], length.out = 8)
#'   x2 <- seq(Xbounds[2, 1], Xbounds[2, 2], length.out = 8)
#'   Xnew <- expand.grid(x1 = x1, x2 = x2)
#'   predict(model2, Xnew = Xnew, l = 12)
#' }
#' @importFrom stats optimise
#' @export
predict.TwinBKP <- function(
    object,
    Xnew,
    l         = 50L,
    v         = NULL,
    CI_level  = 0.95,
    threshold = 0.5,
    ...
) {

  # ======================================================
  # 第零步：从 object 中提取训练时保存的所有信息
  # ======================================================

  Xnorm        <- object$Xnorm         # 全训练集归一化坐标 (n x d)
  y            <- object$y             # 全训练集响应 (n x 1)
  m_train      <- object$m             # 全训练集试验次数 (n x 1)
  Xbounds      <- object$Xbounds       # 输入空间边界 (d x 2)
  global_idx   <- object$global_idx    # 全局点行索引
  Xnorm_global <- object$Xnorm_global  # 全局点归一化坐标
  y_global     <- object$y_global
  m_global     <- object$m_global
  theta_global <- object$theta_global  # 全局核超参数
  kernel       <- object$kernel        # 核函数类型
  isotropic    <- object$isotropic
  prior        <- object$prior
  r0           <- object$r0
  p0           <- object$p0
  loss_type    <- object$loss
  g            <- object$g             # 全局点数量

  n <- nrow(Xnorm)   # 总训练样本量
  d <- ncol(Xnorm)   # 输入维度

  # ======================================================
  # 第一步：输入检查与归一化
  # ======================================================

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

  # 归一化 Xnew 到 [0,1]^d
  Xnew_norm <- sweep(Xnew, 2, Xbounds[, 1], "-")
  Xnew_norm <- sweep(Xnew_norm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

  n_new <- nrow(Xnew_norm)  # 预测点数量

  # ======================================================
  # 第二步：l 和 v 参数检查
  # ======================================================

  if (!is.numeric(l) || length(l) != 1 || l <= 0) stop("'l' must be a positive integer.")
  l <- as.integer(l)
  l_eff <- min(l, n)  # 不超过训练集大小

  if (is.null(v)) v <- 2L * g  # 默认 v = 2g
  if (!is.numeric(v) || length(v) != 1 || v <= 0) stop("'v' must be a positive integer.")
  v <- as.integer(v)

  # ======================================================
  # 第三步：构建 kd-tree（用全训练集归一化坐标）
  # ======================================================

  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package 'FNN' is required for kd-tree search. Install with install.packages('FNN').")
  }

  # 对所有预测点一次性做 kNN 查询，返回 n_new x l_eff 的索引矩阵
  knn_result <- FNN::get.knnx(
    data  = Xnorm,      # 训练集作为搜索库
    query = Xnew_norm,  # 预测点作为查询
    k     = l_eff,
    algorithm = "kd_tree"
  )
  local_idx_mat  <- knn_result$nn.index   # (n_new x l_eff) 局部点索引
  local_dist_mat <- knn_result$nn.dist    # (n_new x l_eff) 到各邻居的距离

  # ======================================================
  # 第四步：按公式(14)计算每个预测点的 Wendland 核超参数 theta_l
  #
  # 公式(14)：theta_l = min{ rho : X_n ⊆ ∪ B_rho(x_i), x_i ∈ X_g }
  # 即：覆盖半径 = max_{x ∈ X_n} min_{x_i ∈ X_g} ||x - x_i||_2
  # 这是一个与预测点无关的全局量，可以预计算一次
  # ======================================================

  # 计算覆盖半径：对全训练集中每个点，找它到全局点集的最近距离
  knn_global <- FNN::get.knnx(
    data  = Xnorm_global,  # 全局点作为搜索库
    query = Xnorm,         # 全训练集作为查询
    k     = 1L,
    algorithm = "kd_tree"
  )
  # theta_l = 覆盖半径 = 最大的"最近全局点距离"
  theta_l <- max(knn_global$nn.dist)

  # q = floor(d/2) + 3（Wendland 核的阶数参数）
  q_wend <- floor(d / 2) + 3

  # ======================================================
  # 第五步：构建 Validation set（从非全局点中随机选 v 个）
  # ======================================================

  # 非全局点索引
  non_global_idx <- setdiff(seq_len(n), global_idx)

  # 若非全局点数量不足，补充全局点
  v_eff <- min(v, length(non_global_idx))
  if (v_eff < v) {
    warning(sprintf(
      "Only %d non-global points available; validation set size reduced from %d to %d.",
      length(non_global_idx), v, v_eff
    ))
  }

  # 随机选取 v_eff 个 validation 点
  val_idx   <- sample(non_global_idx, size = v_eff, replace = FALSE)
  Xnorm_val <- Xnorm[val_idx, , drop = FALSE]
  y_val     <- y[val_idx, , drop = FALSE]
  m_val     <- m_train[val_idx, , drop = FALSE]

  # ======================================================
  # 第六步：预计算全局核矩阵（在 validation 点上，只需计算一次）
  #
  # K_g(val, global)：validation 点 × 全局点 的核矩阵
  # ======================================================

  K_g_val <- kernel_matrix(
    X = Xnorm_val,
    theta = theta_global,
    kernel = kernel,
    isotropic = isotropic
  )

  # ======================================================
  # 第七步：对每个预测点，优化混合比例 lambda
  # ======================================================

  # Wendland 核函数（紧支撑，保证局部核和全局核可识别）
  # L(xa, xb) = (q * r/theta + 1) * max(0, 1 - r/theta)^q
  # 其中 r = ||xa - xb||_2，q = floor(d/2) + 3
  .wendland_kernel <- function(X1, X2, theta) {
    n1 <- nrow(X1); n2 <- nrow(X2)
    K  <- matrix(0, n1, n2)
    for (i in seq_len(n1)) {
      r <- sqrt(rowSums((X2 - matrix(X1[i, ], nrow = n2, ncol = ncol(X1), byrow = TRUE))^2))
      u <- r / theta
      K[i, ] <- (q_wend * u + 1) * pmax(0, 1 - u)^q_wend
    }
    K
  }

  # 混合核 loss（在 validation set 上评估）
  # 混合核矩阵 = lambda * K_g + (1-lambda) * K_l
  .mixed_loss <- function(lambda, K_g, K_l, y_v, m_v, alpha0_v, beta0_v) {
    lambda  <- pmin(pmax(lambda, 0), 1)  # 限制在 [0,1]
    K_mix <- lambda * K_g_val + (1 - lambda) * K_l_val

    if (loss_type == "brier") {
      loss_fun_brier_bkp_rcpp(K_mix, as.numeric(y_v), as.numeric(m_v),
                               as.numeric(alpha0_v), as.numeric(beta0_v))
    } else {
      loss_fun_logloss_bkp_rcpp(K_mix, as.numeric(y_v), as.numeric(m_v),
                                 as.numeric(alpha0_v), as.numeric(beta0_v))
    }
  }

  # ======================================================
  # 第八步：对每个预测点逐一预测
  # ======================================================

  # 预分配输出容器
  pred_mean     <- numeric(n_new)
  pred_var      <- numeric(n_new)
  pred_lower    <- numeric(n_new)
  pred_upper    <- numeric(n_new)
  pred_lambda   <- numeric(n_new)
  pred_theta_l  <- rep(theta_l, n_new)  # theta_l 是全局量，对所有预测点相同
  local_idx_out <- vector("list", n_new)

  for (i in seq_len(n_new)) {

    # --- 8a. 提取第 i 个预测点的局部点 ---
    loc_idx  <- local_idx_mat[i, ]           # 局部点在训练集中的行索引
    Xnorm_loc <- Xnorm[loc_idx, , drop = FALSE]
    y_loc     <- y[loc_idx, , drop = FALSE]
    m_loc     <- m_train[loc_idx, , drop = FALSE]
    local_idx_out[[i]] <- loc_idx

    # --- 8b. 在 validation 点上计算局部 Wendland 核矩阵 ---
    K_l_val <- .wendland_kernel(
      X1 = Xnorm_val,
      X2 = Xnorm_val,
      theta = theta_l
    )

    # --- 8c. 计算 validation 点的先验参数 ---
    # 这里先验基于 validation 点本身，使用 noninformative 先验简化
    # （若需要 adaptive，可用全局 K 计算）
    prior_val <- get_prior(
      prior = prior, model = "BKP",
      r0 = r0, p0 = p0,
      y = y_val, m = m_val,
      K = K_g_val  # 用全局核做先验（最稳）
    )
    alpha0_val <- prior_val$alpha0
    beta0_val  <- prior_val$beta0

    # --- 8d. 优化 lambda（在 validation set 上） ---
    opt_lambda <- optimise(
      f        = .mixed_loss,
      interval = c(0, 1),
      K_g      = K_g_val,
      K_l      = K_l_val,
      y_v      = y_val,
      m_v      = m_val,
      alpha0_v = alpha0_val,
      beta0_v  = beta0_val,
      maximum  = FALSE
    )
    lambda_i <- opt_lambda$minimum
    pred_lambda[i] <- lambda_i

    # --- 8e. 用混合核对预测点做后验更新 ---
    # 预测点到全局点的核向量（1 x g）
    K_g_star <- kernel_matrix(
      X        = matrix(Xnew_norm[i, ], nrow = 1),
      Xprime   = Xnorm_global,
      theta    = theta_global,
      kernel   = kernel,
      isotropic = isotropic
    )  # (1 x g)

    # 预测点到局部点的 Wendland 核向量（1 x l）
    K_l_star <- .wendland_kernel(
      X1    = matrix(Xnew_norm[i, ], nrow = 1),
      X2    = Xnorm_loc,
      theta = theta_l
    )  # (1 x l)

    # 全局点的先验参数
    prior_g <- get_prior(
      prior = prior, model = "BKP",
      r0 = r0, p0 = p0,
      y = y_global, m = m_global,
      K = kernel_matrix(Xnorm_global, theta = theta_global,
                        kernel = kernel, isotropic = isotropic)
    )

    # 局部点的先验参数
    K_ll <- .wendland_kernel(Xnorm_loc, Xnorm_loc, theta_l)
    prior_l <- get_prior(
      prior = prior, model = "BKP",
      r0 = r0, p0 = p0,
      y = y_loc, m = m_loc,
      K = K_ll
    )

    # 混合后验参数
    # alpha_n = alpha0 + lambda * K_g * y_g + (1-lambda) * K_l * y_l
    alpha_n_i <- lambda_i * (as.numeric(prior_g$alpha0) +
                               as.numeric(K_g_star %*% as.numeric(y_global))) +
      (1 - lambda_i) * (as.numeric(prior_l$alpha0) +
                          as.numeric(K_l_star %*% as.numeric(y_loc)))

    beta_n_i  <- lambda_i * (as.numeric(prior_g$beta0) +
                               as.numeric(K_g_star %*% as.numeric(m_global - y_global))) +
      (1 - lambda_i) * (as.numeric(prior_l$beta0) +
                          as.numeric(K_l_star %*% as.numeric(m_loc - y_loc)))

    # 后验均值和方差（Beta 分布）
    eps <- 1e-10
    ab_sum <- max(alpha_n_i + beta_n_i, eps)

    pred_mean[i] <- alpha_n_i / ab_sum
    pred_mean[i] <- min(max(pred_mean[i], eps), 1 - eps)
    pred_var[i]  <- pred_mean[i] * (1 - pred_mean[i]) / (ab_sum + 1)

    # 置信区间（Beta 分位数）
    pred_lower[i] <- suppressWarnings(qbeta((1 - CI_level) / 2, alpha_n_i, beta_n_i))
    pred_upper[i] <- suppressWarnings(qbeta((1 + CI_level) / 2, alpha_n_i, beta_n_i))
  }

  # ======================================================
  # 第九步：组装并返回结果
  # ======================================================

  result <- list(
    Xnew      = Xnew,           # 原始预测点坐标
    Xnew_norm = Xnew_norm,      # 归一化预测点坐标
    mean      = pred_mean,      # 后验均值
    variance  = pred_var,       # 后验方差
    lower     = pred_lower,     # 置信区间下界
    upper     = pred_upper,     # 置信区间上界
    lambda    = pred_lambda,    # 每个预测点的混合比例
    theta_l   = pred_theta_l,   # 局部核超参数（覆盖半径）
    theta_global = theta_global,# 全局核超参数
    local_idx = local_idx_out,  # 每个预测点的局部邻居索引
    CI_level  = CI_level,
    l         = l_eff,
    v         = v_eff
  )

  # 分类标签（仅 m=1 的情形）
  if (all(m_train == 1)) {
    result$class     <- ifelse(pred_mean > threshold, 1, 0)
    result$threshold <- threshold
  }

  class(result) <- "predict_TwinBKP"
  return(result)
}