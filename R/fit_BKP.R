#' @name fit.BKP
#'
#' @title Fit a Beta Kernel Process (BKP) Model
#'
#' @description Fits a Beta Kernel Process model to binomial data
#'    using local kernel smoothing of success and failure counts.
#'
#' @param X Covariate matrix of size \eqn{n \times d}.
#' @param y Vector of observed successes.
#' @param m Vector of binomial counts (number of trials), same length as \code{y}.
#' @param Xbounds Optional \eqn{d \times 2} matrix of lower and upper bounds for each input
#'    dimension (used for normalization). Defaults to \eqn{[0,1]^d}.
#' @param prior Prior type: one of \code{"noninformative"}, \code{"fixed"}, or \code{"adaptive"}.
#' @param r0 Global prior precision (only used in fixed or adaptive priors).
#' @param p0 Global prior mean (only used for fixed priors).
#' @param kernel Kernel function: one of \code{"gaussian"}, \code{"matern52"}, or \code{"matern32"}.
#' @param loss Loss function used for hyperparameter tuning: \code{"brier"} (default) or \code{"log_loss"}.
#' @param n_multi_start Number of random initializations for multi-start optimization (default: \code{10 * d}).
#'
#' @details The BKP model fits a smoothed latent success probability function via
#' kernel-weighted updates to the Beta prior. Inputs are normalized to \eqn{[0,1]^d},
#' and kernel hyperparameters are optimized by minimizing a cross-validated loss.
#'
#' Posterior Beta parameters are computed as:
#' \deqn{
#'   \alpha_n(\mathbf{x}) = \alpha_0(\mathbf{x}) + \sum k(\mathbf{x}, \mathbf{x}_i) y_i,\quad
#'   \beta_n(\mathbf{x}) = \beta_0(\mathbf{x}) + \sum k(\mathbf{x}, \mathbf{x}_i)(m_i - y_i)
#' }
#'
#' @return A list of class \code{"BKP"} representing the fitted Beta Kernel Process model.
#' It contains the following components:
#' \describe{
#'   \item{\code{theta_opt}}{Numeric vector of optimized kernel hyperparameters (lengthscales),
#'     obtained by minimizing the specified loss function.}
#'   \item{\code{kernel}}{Character string specifying the kernel function used, one of
#'     \code{"gaussian"}, \code{"matern52"}, or \code{"matern32"}.}
#'   \item{\code{loss}}{Character string indicating the loss function used during optimization,
#'     either \code{"brier"} or \code{"log_loss"}.}
#'   \item{\code{loss_min}}{Numeric value representing the minimum loss achieved over all
#'     optimization runs.}
#'   \item{\code{X}}{Original (unnormalized) input matrix of dimension \code{n × d}.}
#'   \item{\code{Xnorm}}{Normalized input matrix scaled to the unit hypercube \code{[0,1]^d}.}
#'   \item{\code{Xbounds}}{A \code{d × 2} matrix specifying the lower and upper bounds used
#'     for normalizing each input dimension.}
#'   \item{\code{y}}{Numeric vector of observed binomial successes, of length \code{n}.}
#'   \item{\code{m}}{Numeric vector of binomial trial counts, of length \code{n}.}
#'   \item{\code{prior}}{Character string indicating the prior type used:
#'     \code{"noninformative"}, \code{"fixed"}, or \code{"adaptive"}.}
#'   \item{\code{r0}}{Positive scalar specifying the global precision parameter used
#'     for informative priors.}
#'   \item{\code{p0}}{Prior mean used when \code{prior = "fixed"}.}
#'   \item{\code{alpha0}}{Prior shape parameters \eqn{\alpha_0(\mathbf{x})} for the Beta distribution;
#'     a scalar or length-\code{n} vector.}
#'   \item{\code{beta0}}{Prior shape parameters \eqn{\beta_0(\mathbf{x})} for the Beta distribution;
#'     a scalar or length-\code{n} vector.}
#'   \item{\code{alpha_n}}{Posterior alpha parameters \eqn{\alpha_n(\mathbf{x})}.}
#'   \item{\code{beta_n}}{Posterior beta parameters \eqn{\beta_n(\mathbf{x})}.}
#' }
#'
#' @examples
#' ### 1D
#' set.seed(123)
#' n <- 30
#' Xbounds <- matrix(c(-2,2), nrow=1)
#' x <- tgp::lhs(n = n, rect = Xbounds)
#' m <- sample(100, n, replace = TRUE)
#' true_pi <- (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' y <- rbinom(n, size = m, prob = true_pi)
#' model1 <- fit.BKP(x, y, m, Xbounds=Xbounds)
#' print(model1)
#'
#' ### 2D
#' set.seed(123)
#' n <- 100
#' f <- function(X) {
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
#'   return(f) }
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' x <- tgp::lhs(n = n, rect = Xbounds)
#' m <- sample(100, n, replace = TRUE)
#' true_pi <- pnorm(f(x))
#' y <- rbinom(n, size = m, prob = true_pi)
#' model2 <- fit.BKP(x, y, m)
#' print(model2)
#'
#' @export

fit.BKP <- function(
    X, y, m, Xbounds = NULL,
    prior = c("noninformative", "fixed", "adaptive"), r0 = 2, p0 = 0.5,
    kernel = c("gaussian", "matern52", "matern32"),
    loss = c("brier", "log_loss"),
    n_multi_start = NULL
){
  # ---- Parse and validate arguments ----
  prior <- match.arg(prior)
  kernel <- match.arg(kernel)
  loss <- match.arg(loss)

  # Convert input to proper form
  X <- as.matrix(X)
  y <- matrix(y, ncol = 1)
  m <- matrix(m, ncol = 1)
  d <- ncol(X)
  n <- nrow(X)

  # ---- Validity checks on inputs ----
  if (nrow(y) != n || nrow(m) != n) stop("'y' and 'm' must match nrow(X).")
  if (any(y < 0) || any(m <= 0) || any(y > m)) stop("'y' must be in [0, m] and 'm' > 0.")

  # ---- Normalize input X to [0,1]^d ----
  if (is.null(Xbounds)) Xbounds <- cbind(rep(0, d), rep(1, d))
  if (!all(dim(Xbounds) == c(d, 2))) stop("'Xbounds' must be a d x 2 matrix.")
  Xnorm <- sweep(X, 2, Xbounds[,1], "-")
  Xnorm <- sweep(Xnorm, 2, Xbounds[,2] - Xbounds[,1], "/")

  # ---- Determine initial search space for log10(theta) ----
  # We work in log10(theta) space for numerical stability
  gamma_bounds <- matrix(c((log10(d)-log10(500))/2,       # lower bound
                           (log10(d)+2)/2),               # upper bound
                         ncol = 2, nrow = d, byrow = TRUE)
  if (is.null(n_multi_start)) n_multi_start <- 10 * d
  init_gamma <- lhs(n_multi_start, gamma_bounds)

  # ---- Run multi-start L-BFGS-B optimization to find best kernel parameters ----
  opt_res <- multistart(
    parmat = init_gamma,
    fn     = loss_fun,
    method = "L-BFGS-B",
    lower  = rep(-10, d), # relaxed lower bound
    upper  = rep(10, d),  # relaxed upper bound
    prior = prior, r0 = r0, p0 = p0,
    Xnorm = Xnorm, y = y, m=m,
    loss = loss, kernel = kernel,
    control= list(trace=0))

  # ---- Extract optimal kernel parameters and loss ----
  best_index <- which.min(opt_res$value)
  gamma_opt  <- as.numeric(opt_res[best_index, 1:d])
  theta_opt  <- 10^gamma_opt
  loss_min   <- opt_res$value[best_index]

  # ---- Compute kernel matrix at optimized hyperparameters ----
  K <- kernel_matrix(Xnorm, theta = theta_opt, kernel = kernel)

  # ---- Compute prior parameters (alpha0 and beta0) ----
  prior_par <- get_prior(prior = prior, r0 = r0, p0 = p0, y = y, m = m, K = K)
  alpha0 <- prior_par$alpha0
  beta0  <- prior_par$beta0

  # ---- Compute posterior parameters ----
  alpha_n <- alpha0 + as.vector(K %*% y)
  beta_n  <- beta0 + as.vector(K %*% (m - y))

  # ---- Construct and return the fitted model object ----
  model <- list(
    theta_opt = theta_opt, kernel = kernel,
    loss = loss, loss_min = loss_min,
    X = X, Xnorm = Xnorm, Xbounds = Xbounds, y = y, m = m,
    prior = prior, r0 = r0, p0 = p0, alpha0 = alpha0, beta0 = beta0,
    alpha_n = alpha_n, beta_n = beta_n
  )
  class(model) <- "BKP"
  return(model)
}
