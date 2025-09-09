#' @name fit_BKP
#'
#' @title Fit a Beta Kernel Process (BKP) Model
#'
#' @description Fits a Beta Kernel Process (BKP) model to binary or binomial
#'   response data using local kernel smoothing. The method constructs a
#'   flexible latent probability surface by updating Beta priors with
#'   kernel-weighted observations.
#'
#' @param X A numeric input matrix of size \eqn{n \times d}, where each row
#'   corresponds to a covariate vector.
#' @param y A numeric vector of observed successes (length \code{n}).
#' @param m A numeric vector of total binomial trials (length \code{n}),
#'   corresponding to each \code{y}.
#' @param Xbounds Optional \eqn{d \times 2} matrix specifying the lower and
#'   upper bounds of each input dimension. Used to normalize inputs to
#'   \eqn{[0,1]^d}. If \code{NULL}, inputs are assumed to be pre-normalized, and
#'   default bounds \eqn{[0,1]^d} are applied.
#' @param prior Type of prior: \code{"noninformative"} (default),
#'   \code{"fixed"}, or \code{"adaptive"}.
#' @param r0 Global prior precision (used when \code{prior = "fixed"} or
#'   \code{"adaptive"}).
#' @param p0 Global prior mean (used when \code{prior = "fixed"}). Default is
#'   \code{mean(y/m)}.
#' @param kernel Kernel function for local weighting: \code{"gaussian"}
#'   (default), \code{"matern52"}, or \code{"matern32"}.
#' @param loss Loss function for kernel hyperparameter tuning: \code{"brier"}
#'   (default) or \code{"log_loss"}.
#' @param n_multi_start Number of random initializations for multi-start
#'   optimization. Default is \eqn{10 \times d}.
#' @param theta Optional. A positive scalar or numeric vector of length \code{d}
#'   specifying kernel lengthscale parameters directly. If \code{NULL}
#'   (default), lengthscales are optimized using multi-start L-BFGS-B to
#'   minimize the specified loss.
#'
#' @return A list of class \code{"BKP"} containing the fitted BKP model,
#'   including:
#' \describe{
#'   \item{\code{theta_opt}}{Optimized kernel hyperparameters (lengthscales).}
#'   \item{\code{kernel}}{Kernel function used, as a string.}
#'   \item{\code{loss}}{Loss function used for hyperparameter tuning.}
#'   \item{\code{loss_min}}{Minimum loss achieved during optimization, or
#'     \code{NA} if \code{theta} was user-specified.}
#'   \item{\code{X}}{Original input matrix (\eqn{n \times d}).}
#'   \item{\code{Xnorm}}{Normalized input matrix scaled to \eqn{[0,1]^d}.}
#'   \item{\code{Xbounds}}{Normalization bounds for each input dimension (\eqn{d \times 2}).}
#'   \item{\code{y}}{Observed success counts.}
#'   \item{\code{m}}{Observed binomial trial counts.}
#'   \item{\code{prior}}{Type of prior used.}
#'   \item{\code{r0}}{Prior precision parameter.}
#'   \item{\code{p0}}{Prior mean (for fixed priors).}
#'   \item{\code{alpha0}}{Prior Beta shape parameter \eqn{\alpha_0(\mathbf{x})}.}
#'   \item{\code{beta0}}{Prior Beta shape parameter \eqn{\beta_0(\mathbf{x})}.}
#'   \item{\code{alpha_n}}{Posterior shape parameter \eqn{\alpha_n(\mathbf{x})}.}
#'   \item{\code{beta_n}}{Posterior shape parameter \eqn{\beta_n(\mathbf{x})}.}
#' }
#'
#' @seealso \code{\link{fit_DKP}} for modeling multinomial responses via the
#'   Dirichlet Kernel Process. \code{\link{predict.BKP}},
#'   \code{\link{plot.BKP}}, \code{\link{simulate.BKP}}, and
#'   \code{\link{summary.BKP}} for prediction, visualization, posterior
#'   simulation, and summarization of a fitted BKP model.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv.
#'   https://doi.org/10.48550/arXiv.2508.10447.
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
#' Xbounds <- matrix(c(-2,2), nrow=1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit BKP model
#' model1 <- fit_BKP(X, y, m, Xbounds=Xbounds)
#' print(model1)
#'
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
#' n <- 100
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit BKP model
#' model2 <- fit_BKP(X, y, m, Xbounds=Xbounds)
#' print(model2)
#'
#' @export

fit_BKP <- function(
    X, y, m, Xbounds = NULL,
    prior = c("noninformative", "fixed", "adaptive"), r0 = 2, p0 = mean(y/m),
    kernel = c("gaussian", "matern52", "matern32"),
    loss = c("brier", "log_loss"),
    n_multi_start = NULL, theta = NULL
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

  # ---- prior, kernel, loss ----
  prior  <- match.arg(prior)
  kernel <- match.arg(kernel)
  loss   <- match.arg(loss)

  # ---- Xbounds checks ----
  if (is.null(Xbounds)) {
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

  if (!is.numeric(p0) || length(p0) != 1 || p0 <= 0) {
    stop("'p0' must be a positive scalar.")
  }

  # ---- hyperparameters checks ----
  if (!is.null(n_multi_start)) {
    if (!is.numeric(n_multi_start) || length(n_multi_start) != 1 || n_multi_start <= 0) {
      stop("'n_multi_start' must be a positive integer.")
    }
  }
  if (!is.null(theta)) {
    if (!is.numeric(theta)) stop("'theta' must be numeric.")
    if (length(theta) == 1) {
      theta <- rep(theta, d)
    } else if (length(theta) != d) {
      stop(paste0("'theta' must be either a scalar or a vector of length ", d, "."))
    }
    if (any(theta <= 0)) stop("'theta' must be strictly positive.")
  }

  # ---- Normalize input X to [0,1]^d ----
  Xnorm <- sweep(X, 2, Xbounds[,1], "-")
  Xnorm <- sweep(Xnorm, 2, Xbounds[,2] - Xbounds[,1], "/")

  if (is.null(theta)) {
    # ---- Determine initial search space for log10(theta) ----
    # We work in log10(theta) space for numerical stability
    gamma_bounds <- matrix(c((log10(d) - log10(500))/2,   # lower bound
                             (log10(d) + 2)/2),           # upper bound
                           ncol = 2, nrow = d, byrow = TRUE)
    if (is.null(n_multi_start)) n_multi_start <- 10 * d
    init_gamma <- lhs(n_multi_start, gamma_bounds)

    # ---- Run multi-start L-BFGS-B optimization to find best kernel parameters ----
    opt_res <- multistart(
      parmat = init_gamma,
      fn     = loss_fun,
      method = "L-BFGS-B",
      lower  = rep(-3, d), # relaxed lower bound
      upper  = rep(3, d),  # relaxed upper bound
      prior = prior, r0 = r0, p0 = p0,
      Xnorm = Xnorm, y = y, m=m,
      model = "BKP", loss = loss, kernel = kernel,
      control= list(trace=0))

    # ---- Extract optimal kernel parameters and loss ----
    # opt_res <- opt_res[opt_res$convergence == 0, , drop=FALSE]
    best_index <- which.min(opt_res$value)
    gamma_opt  <- as.numeric(opt_res[best_index, 1:d])
    theta_opt  <- 10^gamma_opt
    loss_min   <- opt_res$value[best_index]
  }else{
    # ---- Use user-provided theta ----
    theta_opt <- theta
    loss_min <- loss_fun(gamma = log10(theta_opt), Xnorm = Xnorm, y = y, m=m,
                         prior = prior, r0 = r0, p0 = p0,
                         model = "BKP", loss = loss, kernel = kernel)
  }

  # ---- Compute kernel matrix at optimized hyperparameters ----
  K <- kernel_matrix(Xnorm, theta = theta_opt, kernel = kernel)

  # ---- Compute prior parameters (alpha0 and beta0) ----
  prior_par <- get_prior(prior = prior, model = "BKP",
                         r0 = r0, p0 = p0, y = y, m = m, K = K)
  alpha0 <- prior_par$alpha0
  beta0  <- prior_par$beta0

  # ---- Compute posterior parameters ----
  alpha_n <- alpha0 + as.vector(K %*% y)
  beta_n  <- beta0 + as.vector(K %*% (m - y))

  # ---- Construct and return the fitted model ----
  BKP_model <- list(
    theta_opt = theta_opt, kernel = kernel,
    loss = loss, loss_min = loss_min,
    X = X, Xnorm = Xnorm, Xbounds = Xbounds, y = y, m = m,
    prior = prior, r0 = r0, p0 = p0, alpha0 = alpha0, beta0 = beta0,
    alpha_n = alpha_n, beta_n = beta_n
  )
  class(BKP_model) <- "BKP"
  return(BKP_model)
}
