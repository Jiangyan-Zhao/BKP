#' @name fit_TwinBKP
#'
#' @title Fit a Twin Beta Kernel Process (TwinBKP) Model
#'
#' @description
#' Fits a TwinBKP model for binary or binomial response data.
#' The workflow follows \code{fit_BKP()}, but kernel hyperparameter tuning
#' is accelerated by first selecting \code{g} representative support points
#' from the full dataset using the \code{Twining} package, then optimizing
#' the kernel lengthscales on these \code{g} points only.
#' The final posterior update still uses all \code{n} observations.
#'
#' @inheritParams fit_BKP
#' @param g Positive integer. Number of global support points selected by
#'   the Twining algorithm for hyperparameter tuning.
#'
#' @return A list of class \code{"TwinBKP"} with the same structure as
#'   \code{"BKP"}, plus:
#' \describe{
#'   \item{\code{g}}{Number of support points used for tuning.}
#'   \item{\code{tune_idx}}{Row indices of the selected support points.}
#' }
#'
#' @seealso \code{\link{fit_BKP}}
#'
#' @examples
#' set.seed(123)
#' true_pi_fun <- function(x) {
#'   (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' }
#' n <- 50
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#' model <- fit_TwinBKP(X, y, m, Xbounds = Xbounds, g = 20)
#' print(model)
#'
#' @export
fit_TwinBKP <- function(
    X, y, m, Xbounds = NULL,
    prior = c("noninformative", "fixed", "adaptive"), r0 = 2, p0 = mean(y / m),
    kernel = c("gaussian", "matern52", "matern32"),
    loss = c("brier", "log_loss"),
    n_multi_start = NULL, theta = NULL,
    isotropic = TRUE,
    g = 20
) {

  # ======================================================
  # 第一步：输入合法性检查
  # ======================================================

  if (missing(X) || missing(y) || missing(m)) {
    stop("Arguments 'X', 'y', and 'm' must be provided.")
  }
  if (!is.matrix(X) && !is.data.frame(X)) stop("'X' must be a numeric matrix or data frame.")
  if (!is.numeric(as.matrix(X))) stop("'X' must contain numeric values only.")
  if (!is.numeric(y)) stop("'y' must be numeric.")
  if (!is.numeric(m)) stop("'m' must be numeric.")

  # 统一转为矩阵列向量格式
  X <- as.matrix(X)
  y <- matrix(y, ncol = 1)
  m <- matrix(m, ncol = 1)

  d <- ncol(X)  # 输入维度
  n <- nrow(X)  # 总样本量

  # 维度一致性检查
  if (nrow(y) != n) stop("'y' must have the same number of rows as 'X'.")
  if (nrow(m) != n) stop("'m' must have the same number of rows as 'X'.")

  # 响应变量合法性
  if (any(y < 0))  stop("'y' must be nonnegative.")
  if (any(m <= 0)) stop("'m' must be strictly positive.")
  if (any(y > m))  stop("Each element of 'y' must be <= corresponding element of 'm'.")
  if (anyNA(X) || anyNA(y) || anyNA(m)) {
    stop("Missing values are not allowed in 'X', 'y', or 'm'.")
  }

  # 解析字符参数
  prior  <- match.arg(prior)
  kernel <- match.arg(kernel)
  loss   <- match.arg(loss)

  # g 合法性检查，限制不超过样本量
  if (!is.numeric(g) || length(g) != 1 || g <= 0) stop("'g' must be a positive integer.")
  g     <- as.integer(g)
  g_eff <- min(g, n)

  # ======================================================
  # 第二步：Xbounds 检查与归一化
  # ======================================================

  if (is.null(Xbounds)) {
    # 未提供 Xbounds 时，检查是否已在 [0,1]
    xmin <- min(X); xmax <- max(X)
    if (xmin < 0 || xmax > 1) {
      warning(sprintf(
        paste0(
          "Input X does not appear to be normalized to [0,1]. ",
          "Current range: [%.3f, %.3f]. ",
          "Please specify Xbounds explicitly."
        ), xmin, xmax
      ))
    }
    Xbounds <- cbind(rep(0, d), rep(1, d))  # 默认 [0,1]^d
  } else {
    if (!is.matrix(Xbounds))   stop("'Xbounds' must be a numeric matrix.")
    if (!is.numeric(Xbounds))  stop("'Xbounds' must contain numeric values.")
    if (!all(dim(Xbounds) == c(d, 2))) {
      stop(paste0("'Xbounds' must be a d x 2 matrix, where d = ", d, "."))
    }
    if (any(Xbounds[, 2] <= Xbounds[, 1])) {
      stop("Each row of 'Xbounds' must satisfy lower < upper.")
    }
  }

  # 将 X 归一化到 [0,1]^d
  Xnorm <- sweep(X, 2, Xbounds[, 1], "-")
  Xnorm <- sweep(Xnorm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

  # ======================================================
  # 第三步：先验与超参数检查
  # ======================================================

  if (!is.numeric(r0) || length(r0) != 1 || r0 <= 0) stop("'r0' must be a positive scalar.")
  if (!is.numeric(p0) || length(p0) != 1 || p0 <= 0 || p0 >= 1) stop("'p0' must be in (0,1).")

  if (!is.null(n_multi_start)) {
    if (!is.numeric(n_multi_start) || length(n_multi_start) != 1 || n_multi_start <= 0) {
      stop("'n_multi_start' must be a positive integer.")
    }
  }

  if (!is.logical(isotropic) || length(isotropic) != 1) {
    stop("'isotropic' must be a single logical value.")
  }

  # 若用户手动指定 theta，校验并扩展
  if (!is.null(theta)) {
    if (!is.numeric(theta)) stop("'theta' must be numeric.")
    if (isotropic && length(theta) != 1) {
      stop("When isotropic=TRUE, 'theta' must be a scalar.")
    }
    if (!isotropic && !(length(theta) == 1 || length(theta) == d)) {
      stop(paste0("When isotropic=FALSE, 'theta' must be scalar or length ", d, "."))
    }
    if (!isotropic && length(theta) == 1) theta <- rep(theta, d)
    if (any(theta <= 0)) stop("'theta' must be strictly positive.")
  }

  # ======================================================
  # 第四步：用 Twining 包选取 g 个全局支撑点
  # ======================================================

  # twin() 返回所选点的行索引（长度为 g_eff）
  tw <- get_twin_indices(Xnorm, g = g_eff, v = 2L * g_eff, runs = 10L, seed = 123L)
  global_idx <- as.integer(tw$gIndices)

  # 提取全局子集数据
  X_global    <- X[global_idx, , drop = FALSE]      # 原始坐标
  Xnorm_global <- Xnorm[global_idx, , drop = FALSE]  # 归一化坐标
  y_global    <- y[global_idx, , drop = FALSE]
  m_global    <- m[global_idx, , drop = FALSE]

  # ======================================================
  # 第五步：在 g 个全局点上优化核超参数（与 fit_BKP 完全相同的逻辑）
  # ======================================================

  if (is.null(theta)) {
    # 超参数维度：各向同性 1 个，各向异性 d 个
    n_theta <- if (isotropic) 1L else d

    # gamma = log10(theta) 的搜索范围
    gamma_bounds <- matrix(
      c((log10(d) - log10(500)) / 2, (log10(d) + 2) / 2),
      ncol = 2, nrow = n_theta, byrow = TRUE
    )

    # 多重启动次数
    if (is.null(n_multi_start)) n_multi_start <- 10L * n_theta

    # LHS 生成多重启动初始点
    init_gamma <- lhs(n_multi_start, gamma_bounds)

    # L-BFGS-B 多重启动优化，目标函数在全局 g 个点上评估
    opt_res <- multistart(
      parmat = init_gamma,
      fn     = loss_fun,
      method = "L-BFGS-B",
      lower  = rep(-3, n_theta),
      upper  = rep(3, n_theta),
      # loss_fun 所需额外参数
      prior = prior, r0 = r0, p0 = p0,
      Xnorm = Xnorm_global,
      y = y_global,
      m = m_global,
      model = "BKP", loss = loss,
      kernel = kernel, isotropic = isotropic,
      control = list(trace = 0)
    )

    # 取所有启动中最优结果
    best_index  <- which.min(opt_res$value)
    gamma_opt   <- as.numeric(opt_res[best_index, 1:n_theta])
    theta_global <- 10^gamma_opt   # 还原到原始尺度
    loss_global  <- opt_res$value[best_index]

  } else {
    # 用户已指定 theta，直接计算损失（不优化）
    theta_global <- theta
    loss_global  <- loss_fun(
      gamma  = log10(theta_global),
      Xnorm  = Xnorm_global,
      y      = y_global,
      m      = m_global,
      prior  = prior, r0 = r0, p0 = p0,
      model  = "BKP", loss = loss,
      kernel = kernel, isotropic = isotropic
    )
  }

  # ======================================================
  # 第六步：组装并返回结果
  # ======================================================

  TwinBKP_model <- list(

    # ---- 全局超参数（核心输出，供 predict_TwinBKP 使用）----
    theta_global = theta_global,  # 全局核超参数
    loss_global  = loss_global,   # 全局点上的最优损失

    # ---- 全局点信息（供 predict_TwinBKP 找局部点时参考）----
    global_idx   = global_idx,    # 全局点在训练集中的行索引
    X_global     = X_global,      # 全局点原始坐标
    Xnorm_global = Xnorm_global,  # 全局点归一化坐标
    y_global     = y_global,      # 全局点响应变量
    m_global     = m_global,      # 全局点试验次数

    # ---- 完整训练数据（predict 阶段需要）----
    X      = X,
    Xnorm  = Xnorm,
    Xbounds = Xbounds,
    y      = y,
    m      = m,

    # ---- 模型配置（predict 阶段继承）----
    kernel    = kernel,
    isotropic = isotropic,
    prior     = prior,
    r0        = r0,
    p0        = p0,
    loss      = loss,
    g         = g_eff
  )

  class(TwinBKP_model) <- "TwinBKP"
  return(TwinBKP_model)
}