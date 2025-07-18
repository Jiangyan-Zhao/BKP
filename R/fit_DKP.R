#' @name fit.DKP
#'
#' @title Fit a Dirichlet Kernel Process (DKP) Model
#'
#' @description Fits a DKP model for multinomial response data by locally
#'   smoothing observed counts to estimate latent Dirichlet parameters.
#'
#' @inheritParams fit.BKP
#' @param Y Matrix of observed multinomial counts, with dimension \eqn{n \times
#'   q}.
#' @param p0 Global prior mean vector (only used when \code{prior = "fixed"}).
#'   Must be of length \eqn{q}.
#'
#' @return A list of class \code{"DKP"} representing the fitted DKP model, with
#'   the following components:
#' \describe{
#'   \item{\code{theta_opt}}{Optimized kernel hyperparameters (lengthscales).}
#'   \item{\code{kernel}}{Kernel function used, as a string.}
#'   \item{\code{loss}}{Loss function used for hyperparameter tuning.}
#'   \item{\code{loss_min}}{Minimum loss value achieved during optimization.}
#'   \item{\code{X}}{Original (unnormalized) input matrix of size \code{n × d}.}
#'   \item{\code{Xnorm}}{Normalized input matrix scaled to \eqn{[0,1]^d}.}
#'   \item{\code{Xbounds}}{Matrix specifying normalization bounds for each input dimension.}
#'   \item{\code{Y}}{Observed multinomial counts of size \code{n × q}.}
#'   \item{\code{prior}}{Type of prior used.}
#'   \item{\code{r0}}{Prior precision parameter.}
#'   \item{\code{p0}}{Prior mean (for fixed priors).}
#'   \item{\code{alpha0}}{Prior Dirichlet parameters at each input location (scalar or matrix).}
#'   \item{\code{alpha_n}}{Posterior Dirichlet parameters after kernel smoothing.}
#' }
#'
#' @seealso \code{\link{fit.BKP}} for modeling binomial responses using the Beta
#'   Kernel Process. \code{\link{predict.DKP}}, \code{\link{plot.DKP}},
#'   \code{\link{simulate.DKP}} for making predictions, visualizing results, and
#'   generating simulations from a fitted DKP model. \code{\link{summary.DKP}},
#'   \code{\link{print.DKP}} for inspecting fitted model summaries.
#'
#' @examples
#' #-------------------------- 1D Example ---------------------------
#' set.seed(123)
#'
#' # Define true class probability function (3-class)
#' true_pi_fun <- function(X) {
#'   p <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
#'   return(matrix(c(p/2, p/2, 1 - p), nrow = length(p)))
#' }
#'
#' n <- 30
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit DKP model
#' model1 <- fit.DKP(X, Y, Xbounds = Xbounds)
#' print(model1)
#'
#'
#' #-------------------------- 2D Example ---------------------------
#' set.seed(123)
#'
#' # Define latent function and transform to 3-class probabilities
#' true_pi_fun <- function(X) {
#'   if (is.null(nrow(X))) X <- matrix(X, nrow = 1)
#'   m <- 8.6928; s <- 2.4269
#'   x1 <- 4 * X[,1] - 2
#'   x2 <- 4 * X[,2] - 2
#'   a <- 1 + (x1 + x2 + 1)^2 *
#'     (19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2)
#'   b <- 30 + (2*x1 - 3*x2)^2 *
#'     (18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2)
#'   f <- (log(a * b) - m) / s
#'   p <- pnorm(f)
#'   return(matrix(c(p/2, p/2, 1 - p), nrow = length(p)))
#' }
#'
#' n <- 100
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit DKP model
#' model2 <- fit.DKP(X, Y, Xbounds = Xbounds)
#' print(model2)
#'
#' @export

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
