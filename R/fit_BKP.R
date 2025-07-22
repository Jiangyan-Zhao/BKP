#' @name fit.BKP
#'
#' @title Fit a Beta Kernel Process (BKP) Model
#'
#' @description Fits a BKP model to binomial or binary response data via local
#'   kernel smoothing. The model constructs a flexible latent probability
#'   surface by updating Beta priors using kernel-weighted observations.
#'
#' @param X A numeric input matrix of size \eqn{n \times d}, where each row
#'   represents a covariate vector.
#' @param y A numeric vector of observed successes (length \code{n}).
#' @param m A numeric vector of total binomial trials (length \code{n}),
#'   corresponding to each \code{y}.
#' @param Xbounds Optional \eqn{d \times 2} matrix specifying the lower and
#'   upper bounds of each input dimension. Used to normalize inputs to
#'   \eqn{[0,1]^d}. If \code{Xbounds} is \code{NULL}, the input is assumed to
#'   have already been normalized, and the default bounds are set to \eqn{[0,1]^d}.
#' @param prior Type of prior to use. One of \code{"noninformative"},
#'   \code{"fixed"}, or \code{"adaptive"}.
#' @param r0 Global prior precision (only used when \code{prior = "fixed"} or
#'   \code{"adaptive"}).
#' @param p0 Global prior mean (only used when \code{prior = "fixed"}).
#' @param kernel Kernel function for local weighting. Choose from
#'   \code{"gaussian"}, \code{"matern52"}, or \code{"matern32"}.
#' @param loss Loss function for kernel hyperparameter tuning. One of
#'   \code{"brier"} (default) or \code{"log_loss"}.
#' @param n_multi_start Number of random initializations for multi-start
#'   optimization. Default is \code{10 × d}.
#'
#' @return A list of class \code{"BKP"} containing the fitted BKP model, with
#'   the following elements:
#' \describe{
#'   \item{\code{theta_opt}}{Optimized kernel hyperparameters (lengthscales).}
#'   \item{\code{kernel}}{Kernel function used, as a string.}
#'   \item{\code{loss}}{Loss function used for hyperparameter tuning.}
#'   \item{\code{loss_min}}{Minimum loss value achieved during optimization.}
#'   \item{\code{X}}{Original (unnormalized) input matrix of size \code{n × d}.}
#'   \item{\code{Xnorm}}{Normalized input matrix scaled to \eqn{[0,1]^d}.}
#'   \item{\code{Xbounds}}{Matrix specifying normalization bounds for each input dimension.}
#'   \item{\code{y}}{Observed success counts.}
#'   \item{\code{m}}{Observed binomial trial counts.}
#'   \item{\code{prior}}{Type of prior used.}
#'   \item{\code{r0}}{Prior precision parameter.}
#'   \item{\code{p0}}{Prior mean (for fixed priors).}
#'   \item{\code{alpha0}}{Prior shape parameter \eqn{\alpha_0(\mathbf{x})}, either a scalar or vector.}
#'   \item{\code{beta0}}{Prior shape parameter \eqn{\beta_0(\mathbf{x})}, either a scalar or vector.}
#'   \item{\code{alpha_n}}{Posterior shape parameter \eqn{\alpha_n(\mathbf{x})}.}
#'   \item{\code{beta_n}}{Posterior shape parameter \eqn{\beta_n(\mathbf{x})}.}
#' }
#'
#' @seealso \code{\link{fit.DKP}} for modeling multinomial responses using the
#'   Dirichlet Kernel Process. \code{\link{predict.BKP}},
#'   \code{\link{plot.BKP}}, \code{\link{simulate.BKP}} for making predictions,
#'   visualizing results, and generating simulations from a fitted BKP model.
#'   \code{\link{summary.BKP}}, \code{\link{print.BKP}} for inspecting model
#'   details.
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
#' model1 <- fit.BKP(X, y, m, Xbounds=Xbounds)
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
#' model2 <- fit.BKP(X, y, m, Xbounds=Xbounds)
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
