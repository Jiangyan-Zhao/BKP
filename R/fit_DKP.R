#' @name fit.DKP
#'
#' @title Fit a Dirichlet Kernel Process (DKP) Model for Multinomial Observations
#'
#' @description
#' Fits a Dirichlet Kernel Process model for multinomiall data, estimating latent functions
#' representing Dirichlet distribution parameters using kernel smoothing.
#'
#' @param data Optional data frame containing covariates \code{X} and observations \code{Y}
#'   in the first \eqn{d} and \eqn{d+q} columns, respectively.
#' @param X Covariate matrix of size \eqn{n \times d}. Ignored if \code{data} is provided.
#' @param Y Matrix of observed successes of size \eqn{n \times q}.
#' @param Xbounds Optional \eqn{d \times 2} matrix defining lower and upper bounds for each input dimension
#'   (used for normalization). Defaults to \code{[0,1]^d}.
#' @param prior Type of prior used: "noninformative", "fixed", or "adaptive".
#' @param r0 Global precision parameter used in prior specification.
#' @param p0 Global prior mean used when prior = "fixed". Must be a vector of length \eqn{q}.
#' @param kernel Character string specifying the kernel function to use: one of \code{"gaussian"},
#'   \code{"matern52"}, or \code{"matern32"}.
#' @param loss Loss function to be minimized during hyperparameter tuning. Choose between \code{"brier"}
#'   (default) and \code{"NLML"} (negative log marginal likelihood).
#' @param n_multi_start Number of initial points for multi-start optimization (default = 10*d).
#'
#' @details
#' The function fits the DKP model and stores everything necessary for prediction.
#' The kernel hyper-parameters optimization uses multi-start gradient-based method.
#' Inputs are first normalized to the \eqn{[0,1]^d} space defined by \code{Xbounds}.
#' The kernel matrix is computed with optimized parameters, and posterior Dirichlet parameters
#' \eqn{\alpha_n} is computed by smoothing the observed successes \eqn{y} with the kernel matrix:
#' \deqn{
#' \alpha_n = \alpha_0 + K \cdot Y,
#' }
#' where \eqn{K} is the kernel matrix evaluated at normalized inputs.
#'
#' @return A list of class \code{"DKP"} containing the fitted model, with the following components:
#' \describe{
#'   \item{\code{bestTheta}}{A numeric vector of optimal kernel hyperparameters (lengthscale), obtained by minimizing the specified loss function.}
#'   \item{\code{kernel}}{Character string indicating the kernel type used ("gaussian", "matern52", or "matern32").}
#'   \item{\code{loss}}{Character string specifying the loss function used ("brier" or "log_loss").}
#'   \item{\code{minLoss}}{The minimum loss value achieved during hyperparameter optimization.}
#'   \item{\code{X}}{The original input matrix (unnormalized), with dimensions \code{n x d}.}
#'   \item{\code{Xnorm}}{The normalized input matrix scaled to \code{[0,1]^d}.}
#'   \item{\code{Xbounds}}{A \code{d x 2} matrix specifying the lower and upper bounds used for normalizing each dimension of \code{X}.}
#'   \item{\code{Y}}{A numeric matrix of observed successes (length \code{n x q}).}
#'   \item{\code{prior}}{Character string indicating the prior type used ("noninformative", "fixed", or "adaptive").}
#'   \item{\code{r0}}{Global prior precision parameter (only used for fixed or adaptive priors).}
#'   \item{\code{p0}}{Global prior mean parameter (only used for fixed priors).}
#'   \item{\code{alpha0}}{The prior alpha parameters for each input location. A vector or matrix with \code{n} rows.}
#'   \item{\code{alpha_n}}{Posterior alpha parameters computed as \code{alpha0 + K \%*\% Y}.}
#' }
#'
#' @author Jiangyan Zhao, Kunhai Qing, Jin Xu
#' @examples
#' ### 1D
#' set.seed(123)
#' n <- 30
#' Xbounds <- matrix(c(-2,2), nrow=1)
#' x <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' true_pi <- matrix(c(true_pi/2,true_pi/2,1-true_pi),nrow = n, byrow = FALSE)
#' m <- sample(100, n, replace = TRUE)
#' Y <- matrix(0, nrow = n, ncol = 3)
#' for (i in 1:n) {
#'   Y[i, ] <- rmultinom(n=1, size=m[i], prob=true_pi[i, ])
#' }
#' DKPmodel <- fit.DKP(p=3,X=x,Y=Y,Xbounds = Xbounds,prior = "noninformative",kernel = "gaussian",loss = "brier")
#' print(DKPmodel)
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
#' true_pi <- pnorm(f(x))
#' true_pi <- matrix(c(true_pi/2,true_pi/2,1-true_pi),nrow = n, byrow = FALSE)
#' m <- sample(100, n, replace = TRUE)
#' Y <- matrix(0, nrow = n, ncol = 3)
#' for (i in 1:n) {
#'   Y[i, ] <- rmultinom(n=1, size=m[i], prob=true_pi[i, ])
#' }
#' DKPmodel <- fit.DKP(p=3,X=x,Y=Y,Xbounds = Xbounds,prior = "noninformative",kernel = "gaussian",loss = "brier")
#' print(DKPmodel)
#'
#' @export
#' @importFrom optimx multistart

fit.DKP <- function(
    X, Y, Xbounds = NULL,
    prior = c("noninformative", "fixed", "adaptive"), r0 = 2, p0 = NULL,
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
  Y <- as.matrix(Y)
  d <- ncol(X)
  q <- ncol(Y)
  n <- nrow(X)

  # ---- Validity checks on inputs ----
  if (nrow(Y) != n) stop("Number of rows in 'Y' must match number of rows in 'X'.")
  if (any(Y < 0)) stop("'Y' must be in non-negtive.")

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
    fn     = loss_fun_dkp,
    method = "L-BFGS-B",
    lower  = rep(-10, d), # relaxed lower bound
    upper  = rep(10, d),  # relaxed upper bound
    prior = prior, r0 = r0, p0 = p0,
    Xnorm = Xnorm, Y = Y,
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
  alpha0 <- get_prior_dkp(prior = prior, r0 = r0, p0 = p0, Y = Y, K = K)

  # ---- Compute posterior parameters ----
  alpha_n <- alpha0 + as.matrix(K %*% Y)

  # ---- Construct and return the fitted model object ----
  model <- list(
    theta_opt = theta_opt, kernel = kernel,
    loss = loss, loss_min = loss_min,
    X = X, Xnorm = Xnorm, Xbounds = Xbounds, Y = Y,
    prior = prior, r0 = r0, p0 = p0,
    alpha0 = alpha0, alpha_n = alpha_n
  )
  class(model) <- "DKP"
  return(model)
}
