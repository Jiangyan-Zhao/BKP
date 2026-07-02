#' @name fit_BKP
#'
#' @title Fit a Beta Kernel Process (BKP) Model
#'
#' @description Fit a Beta Kernel Process (BKP) model for binary or binomial
#'   response data. The model estimates a covariate-dependent success
#'   probability surface by combining kernel-based smoothing with conjugate Beta
#'   posterior updates.
#'
#' @param X A numeric input matrix or data frame of size \eqn{n \times d}, where
#'   each row is an input or covariate vector.
#' @param y A numeric vector of observed success counts of length \code{n}.
#' @param m A numeric vector of total binomial trial counts of length \code{n}.
#'   Each element of \code{y} must be between 0 and the corresponding element of
#'   \code{m}.
#' @param Xbounds Optional \eqn{d \times 2} numeric matrix giving the lower and
#'   upper bounds of each input dimension. When supplied, \code{X} is normalized
#'   to \eqn{[0,1]^d} before fitting. If \code{NULL}, \code{X} is assumed to be
#'   already normalized to \eqn{[0,1]^d}.
#' @param prior Type of prior specification. Options are
#'   \code{"noninformative"} (default), \code{"fixed"}, and
#'   \code{"adaptive"}.
#' @param r0 Positive scalar prior precision used when \code{prior = "fixed"} or
#'   \code{prior = "adaptive"}. Default is \code{2}.
#' @param p0 Prior mean used when \code{prior = "fixed"}. It must be a scalar in
#'   \eqn{(0,1)}. The default is the empirical mean \code{mean(y / m)}.
#' @param kernel Kernel function used for local weighting. Options are
#'   \code{"gaussian"} (default), \code{"matern52"}, \code{"matern32"}, and
#'   \code{"wendland"}.
#' @param loss Leave-one-out loss used for kernel hyperparameter tuning. Options
#'   are \code{"brier"} (default) and \code{"log_loss"}.
#' @param n_multi_start Number of initial points used in multi-start
#'   derivative-free optimization of the kernel lengthscale parameters. If
#'   \code{NULL}, the default is \eqn{10p}, where \eqn{p = 1} for isotropic
#'   kernels and \eqn{p = d} for anisotropic kernels.
#' @param theta Optional positive kernel lengthscale parameter. When
#'   \code{isotropic = TRUE}, this must be a scalar. When
#'   \code{isotropic = FALSE}, this can be either a scalar, which is broadcast
#'   to all dimensions, or a numeric vector of length \code{d}. If \code{NULL},
#'   the lengthscale parameters are selected by minimizing the specified
#'   leave-one-out loss.
#' @param isotropic Logical. If \code{TRUE} (default), use a single shared
#'   lengthscale across input dimensions. If \code{FALSE}, use separate
#'   per-dimension lengthscales.
#' @param n_threads Number of OpenMP threads used for multi-start optimization.
#'   Default is \code{1}. This argument only affects fitting when
#'   \code{theta = NULL}.
#' @param ess Effective-sample-size calibration for the kernel-weighted data
#'   contribution. Use \code{"none"} (default) for the standard BKP posterior
#'   update. Use \code{"shepard"} to rescale the kernel-weighted data
#'   contribution so that its effective trial size is
#'   \eqn{\rho(\mathbf{x}) m_S(\mathbf{x})}, where \eqn{m_S(\mathbf{x})} is a
#'   Shepard interpolation of the observed trial sizes on the normalized input
#'   scale and \eqn{\rho(\mathbf{x}) = \max_i k(\mathbf{x}, \mathbf{x}_i)}. This
#'   calibration preserves the kernel-weighted empirical proportion and changes
#'   only the data precision, not the prior parameters.
#'
#' @details Inputs are normalized to \eqn{[0,1]^d} using \code{Xbounds}. For a
#'   location \eqn{\mathbf{x}}, BKP computes kernel weights
#'   \eqn{k(\mathbf{x}, \mathbf{x}_i)} between \eqn{\mathbf{x}} and the training
#'   inputs. These weights are used to update a local Beta prior with
#'   kernel-weighted binomial counts:
#'   \deqn{
#'     \alpha_n(\mathbf{x}) =
#'     \alpha_0(\mathbf{x}) +
#'     \sum_i k(\mathbf{x}, \mathbf{x}_i) y_i,
#'   }
#'   \deqn{
#'     \beta_n(\mathbf{x}) =
#'     \beta_0(\mathbf{x}) +
#'     \sum_i k(\mathbf{x}, \mathbf{x}_i) (m_i - y_i).
#'   }
#'
#'   If \code{theta = NULL}, the kernel lengthscale parameters are selected by
#'   leave-one-out cross-validation using the specified \code{loss}. Optimization
#'   is performed over log-transformed lengthscales using a multi-start
#'   derivative-free search. If \code{theta} is supplied, optimization is skipped
#'   and the model is fitted using the supplied lengthscale parameter.
#'
#'   If \code{ess = "shepard"}, only the kernel-weighted data contribution is
#'   rescaled to match the local effective-sample-size target; the prior
#'   parameters are not rescaled. Shepard calibration requires unique training
#'   input locations on the normalized input scale.
#'
#'   The returned object stores posterior parameters evaluated at the training
#'   inputs. Posterior inference at new inputs is performed by
#'   \code{\link{predict.BKP}}.
#'
#' @return A list of class \code{"BKP"} containing the fitted model, with the
#'   following components:
#' \describe{
#'   \item{\code{theta_opt}}{Optimized or user-specified kernel lengthscale
#'     parameter(s).}
#'   \item{\code{kernel}}{Kernel function used.}
#'   \item{\code{isotropic}}{Logical flag indicating whether a shared
#'     lengthscale or per-dimension lengthscales were used.}
#'   \item{\code{loss}}{Loss function used for hyperparameter tuning.}
#'   \item{\code{loss_min}}{Loss value at the selected or user-specified
#'     lengthscale parameter(s).}
#'
#'   \item{\code{X}}{Original input matrix.}
#'   \item{\code{Xnorm}}{Input matrix normalized to \eqn{[0,1]^d}.}
#'   \item{\code{Xbounds}}{Normalization bounds for each input dimension.}
#'   \item{\code{y}}{Observed success counts, stored as a one-column matrix.}
#'   \item{\code{m}}{Observed binomial trial counts, stored as a one-column
#'     matrix.}
#'
#'   \item{\code{prior}}{Prior specification used.}
#'   \item{\code{r0}}{Prior precision parameter.}
#'   \item{\code{p0}}{Prior mean used when \code{prior = "fixed"}.}
#'   \item{\code{alpha0}}{Prior Beta \eqn{\alpha} shape parameters evaluated at
#'     the training inputs.}
#'   \item{\code{beta0}}{Prior Beta \eqn{\beta} shape parameters evaluated at
#'     the training inputs.}
#'
#'   \item{\code{alpha_n}}{Posterior Beta \eqn{\alpha} shape parameters
#'     evaluated at the training inputs.}
#'   \item{\code{beta_n}}{Posterior Beta \eqn{\beta} shape parameters evaluated
#'     at the training inputs.}
#'
#'   \item{\code{ess}}{Effective-sample-size calibration method used.}
#'   \item{\code{ess_info}}{ESS calibration diagnostics, or \code{NULL} when
#'     \code{ess = "none"}.}
#' }
#'
#' @seealso \code{\link{fit_DKP}} for Dirichlet Kernel Process modeling of
#'   categorical or multinomial responses; \code{\link{fit_TwinBKP}} for the
#'   scalable global-local TwinBKP approximation; \code{\link{predict.BKP}},
#'   \code{\link{plot.BKP}}, \code{\link{simulate.BKP}}, and
#'   \code{\link{summary.BKP}} for downstream methods.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv. <doi:10.48550/arXiv.2508.10447>.
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
#' n <- 30
#' Xbounds <- matrix(c(-2,2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit BKP model
#' # A fixed theta is used here only to keep the example fast and reproducible.
#' # In practice, omit theta to select it by leave-one-out cross-validation.
#' model1 <- fit_BKP(X, y, m, Xbounds = Xbounds, theta = 0.04)
#' print(model1)
#'
#'
#' #-------------------------- 2D Example ---------------------------
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
#' # A fixed theta is used here only to keep the example fast and reproducible.
#' # In practice, omit theta to select it by leave-one-out cross-validation.
#' model2 <- fit_BKP(X, y, m, Xbounds = Xbounds, theta = 0.08)
#' print(model2)
#'
#' @export

fit_BKP <- function(
    X, y, m, Xbounds = NULL,
    prior = c("noninformative", "fixed", "adaptive"), r0 = 2, p0 = mean(y/m),
    kernel = c("gaussian", "matern52", "matern32", "wendland"),
    loss = c("brier", "log_loss"),
    n_multi_start = NULL, theta = NULL,
    isotropic = TRUE, n_threads = 1,
    ess = c("none", "shepard")
){
  # ---- Argument checking ----
  if (missing(X) || missing(y) || missing(m)) {
    stop("Arguments 'X', 'y', and 'm' must be provided.")
  }
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("'X' must be a numeric matrix or data frame.")
  }
  if (!is.numeric(as.matrix(X))) {
    stop("'X' must contain numeric values only.")
  }
  if (!is.numeric(y)) stop("'y' must be numeric.")
  if (!is.numeric(m)) stop("'m' must be numeric.")

  X <- as.matrix(X)
  y <- matrix(y, ncol = 1)
  m <- matrix(m, ncol = 1)

  d <- ncol(X)
  n <- nrow(X)

  if (nrow(y) != n) stop("'y' must have the same number of rows as 'X'.")
  if (nrow(m) != n) stop("'m' must have the same number of rows as 'X'.")
  if (any(y < 0)) stop("'y' must be nonnegative.")
  if (any(m <= 0)) stop("'m' must be strictly positive.")
  if (any(y > m)) stop("Each element of 'y' must be less than or equal to corresponding element of 'm'.")
  if (anyNA(X) || anyNA(y) || anyNA(m)) stop("Missing values are not allowed in 'X', 'y', or 'm'.")
  if (any(!is.finite(X)) || any(!is.finite(y)) || any(!is.finite(m))) {
    stop("'X', 'y', and 'm' must contain only finite values.")
  }

  # ---- prior, kernel, loss ----
  prior  <- match.arg(prior)
  kernel <- match.arg(kernel)
  loss   <- match.arg(loss)
  ess    <- match.arg(ess)

  # ---- Xbounds checks ----
  if (is.null(Xbounds)) {
    # Check if X already seems normalized
    xmin <- min(X)
    xmax <- max(X)

    if (xmin < 0 || xmax > 1) {
      warning(
        sprintf(
          paste0(
            "Input X does not appear to be normalized to [0,1]. ",
            "Current range: [%.3f, %.3f].\n",
            "Please normalize X or specify Xbounds explicitly; ",
            "otherwise the model may produce incorrect results."
          ),
          xmin, xmax
        )
      )
    }
    # Default bounds: assume X already in [0,1]^d
    Xbounds <- cbind(rep(0, d), rep(1, d))
  } else {
    if (!is.matrix(Xbounds)) stop("'Xbounds' must be a numeric matrix.")
    if (!is.numeric(Xbounds)) stop("'Xbounds' must contain numeric values.")
    if (!all(dim(Xbounds) == c(d, 2))) {
      stop(paste0("'Xbounds' must be a matrix with dimensions d x 2, where d = ", d, "."))
    }
    if (any(Xbounds[,2] <= Xbounds[,1])) {
      stop("Each row of 'Xbounds' must satisfy lower < upper.")
    }
  }

  # ---- prior parameters checks ----
  if (!is.numeric(r0) || length(r0) != 1 || r0 <= 0) {
    stop("'r0' must be a positive scalar.")
  }

  if (prior == "fixed") {
    if (!is.numeric(p0) || length(p0) != 1 ||
        is.na(p0) || !is.finite(p0) || p0 <= 0 || p0 >= 1) {
      stop("For fixed prior in BKP, 'p0' must be a scalar in (0, 1).")
    }
  }

  # ---- hyperparameters checks ----
  if (!is.null(n_multi_start)) {
    if (!is.numeric(n_multi_start) || length(n_multi_start) != 1  ||
        is.na(n_multi_start) || !is.finite(n_multi_start) || n_multi_start <= 0) {
      stop("'n_multi_start' must be a positive integer.")
    }
  }

  if (!is.numeric(n_threads) || length(n_threads) != 1 ||
      is.na(n_threads) || !is.finite(n_threads) || n_threads <= 0) {
    stop("'n_threads' must be a positive integer.")
  }
  n_threads <- as.integer(n_threads)

  if (!is.null(theta)) {
    if (!is.numeric(theta)) stop("'theta' must be numeric.")
    if (!is.logical(isotropic) || length(isotropic) != 1) {
      stop("'isotropic' must be a single logical value.")
    }
    if (isotropic) {
      if (length(theta) != 1) {
        stop("When isotropic=TRUE, 'theta' must be a scalar.")
      }
    } else if (!(length(theta) == 1 || length(theta) == d)) {
      stop(paste0("When isotropic=FALSE, 'theta' must be either a scalar or a vector of length ", d, "."))
    }
    if (!isotropic && length(theta) == 1) theta <- rep(theta, d)
    if (any(theta <= 0)) stop("'theta' must be strictly positive.")
  } else {
    if (!is.logical(isotropic) || length(isotropic) != 1) {
      stop("'isotropic' must be a single logical value.")
    }
  }

  # ---- Normalize input X to [0,1]^d ----
  Xnorm <- sweep(X, 2, Xbounds[,1], "-")
  Xnorm <- sweep(Xnorm, 2, Xbounds[,2] - Xbounds[,1], "/")

  if (identical(ess, "shepard")) {
    bkp_check_unique_locations(Xnorm)
  }

  if (is.null(theta)) {
    # ---- Determine the number of optimization variables ----
    n_theta <- ifelse(isTRUE(isotropic), 1L, d)

    # ---- Number of multi-start initial points ----
    if (is.null(n_multi_start)) {
      n_multi_start <- as.integer(10L * n_theta)
    } else {
      n_multi_start <- as.integer(n_multi_start)
    }

    # ---- Initial search region Omega_0 for log10(theta) ----
    gamma_bounds <- matrix(c((log10(n_theta) - log10(500))/2,   # lower bound
                             (log10(n_theta) + 2)/2),           # upper bound
                           ncol = 2, nrow = n_theta, byrow = TRUE)
    init_gamma <- lhs(n = n_multi_start, rect = gamma_bounds) # tgp::lhs

    # ---- Local optimization region Omega = [-3, 3]^p ----
    lower <- rep(-3, n_theta)
    upper <- rep(3, n_theta)

    max_iter <- min(500L, ceiling(100 * log1p(n_theta)))

    m_shepard_loo <- if (identical(ess, "shepard")) bkp_shepard_m_loo(Xnorm, m, power = 2) else NULL

    opt_cpp <- optimize_bkp_theta_rcpp(
      Xnorm = Xnorm,
      y = y,
      m = m,
      prior = prior,
      r0 = r0,
      p0 = p0,
      loss = loss,
      kernel = kernel,
      isotropic = isotropic,
      init_gamma = init_gamma,
      lower = lower,
      upper = upper,
      max_iter = max_iter,
      n_threads = n_threads,
      ess = ess,
      m_shepard_loo = m_shepard_loo
    )

    gamma_opt <- as.numeric(opt_cpp$gamma_opt)
    theta_opt <- as.numeric(opt_cpp$theta_opt)
    loss_min <- as.numeric(opt_cpp$loss_min)
  } else {
    theta_opt <- theta
    gamma_opt <- log10(theta_opt)
    loss_min <- loss_fun(
      gamma = gamma_opt, Xnorm = Xnorm, y = y, m = m,
      prior = prior, r0 = r0, p0 = p0,
      model = "BKP", loss = loss, kernel = kernel,
      isotropic = isotropic, ess = ess
    )
  }

  # ---- Compute prior and posterior parameters ----
  posterior <- bkp_compute_posterior(
    Xquery_norm = Xnorm, Xtrain_norm = Xnorm, y = y, m = m,
    theta = theta_opt, kernel = kernel, isotropic = isotropic,
    prior = prior, r0 = r0, p0 = p0, ess = ess
  )
  K <- posterior$K
  alpha0 <- posterior$alpha0
  beta0 <- posterior$beta0
  alpha_n <- posterior$alpha_n
  beta_n <- posterior$beta_n
  ess_info <- posterior$ess_info

  # ---- Construct and return the fitted model ----
  BKP_model <- list(
    # Model configuration
    theta_opt = theta_opt,
    kernel = kernel,
    isotropic = isotropic,
    loss = loss,
    loss_min = loss_min,

    # Training data and normalization
    X = X,
    Xnorm = Xnorm,
    Xbounds = Xbounds,
    y = y,
    m = m,

    # Prior specification
    prior = prior,
    r0 = r0,
    p0 = p0,
    alpha0 = alpha0,
    beta0 = beta0,

    # Posterior parameters at training inputs
    alpha_n = alpha_n,
    beta_n = beta_n,

    # Effective-sample-size calibration
    ess = ess,
    ess_info = ess_info
  )
  class(BKP_model) <- "BKP"
  return(BKP_model)
}
