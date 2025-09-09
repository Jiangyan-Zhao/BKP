#' @name fit_DKP
#'
#' @title Fit a Dirichlet Kernel Process (DKP) Model
#'
#' @description Fits a DKP model for categorical or multinomial response data by
#'   locally smoothing observed counts to estimate latent Dirichlet parameters.
#'   The model constructs flexible latent probability surfaces by updating
#'   Dirichlet priors using kernel-weighted observations.
#'
#' @inheritParams fit_BKP
#' @param Y Numeric matrix of observed multinomial counts, with dimension
#'   \eqn{n \times q}.
#' @param p0 Global prior mean vector (used only when \code{prior = "fixed"}).
#'   Defaults to the empirical class proportions \code{colMeans(Y / rowSums(Y))}.
#'   Must have length equal to the number of categories \eqn{q}.
#'
#' @return A list of class \code{"DKP"} representing the fitted DKP model, with
#'   the following components:
#' \describe{
#'   \item{\code{theta_opt}}{Optimized kernel hyperparameters (lengthscales).}
#'   \item{\code{kernel}}{Kernel function used, as a string.}
#'   \item{\code{loss}}{Loss function used for hyperparameter tuning.}
#'   \item{\code{loss_min}}{Minimum loss value achieved during kernel
#'     hyperparameter optimization. Set to \code{NA} if \code{theta} is user-specified.}
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
#' @seealso \code{\link{fit_BKP}} for modeling binomial responses via the Beta
#'   Kernel Process. \code{\link{predict.DKP}}, \code{\link{plot.DKP}},
#'   \code{\link{simulate.DKP}} for prediction, visualization, and posterior
#'   simulation from a fitted DKP model. \code{\link{summary.DKP}},
#'   \code{\link{print.DKP}} for inspecting model summaries.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv.
#'   https://doi.org/10.48550/arXiv.2508.10447.
#'
#' @examples
#' #-------------------------- 1D Example ---------------------------
#' set.seed(123)
#'
#' # Define true class probability function (3-class)
#' true_pi_fun <- function(X) {
#'   p1 <- 1/(1+exp(-3*X))
#'   p2 <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
#'   return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
#' }
#'
#' n <- 30
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(150, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit DKP model
#' model1 <- fit_DKP(X, Y, Xbounds = Xbounds)
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
#'   f <- (log(a*b)- m)/s
#'   p1 <- pnorm(f) # Transform to probability
#'   p2 <- sin(pi * X[,1]) * sin(pi * X[,2])
#'   return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
#' }
#'
#' n <- 100
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(150, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit DKP model
#' model2 <- fit_DKP(X, Y, Xbounds = Xbounds)
#' print(model2)
#'
#' @export

fit_DKP <- function(
    X, Y, Xbounds = NULL,
    prior = c("noninformative", "fixed", "adaptive"), r0 = 2, p0 = colMeans(Y / rowSums(Y)),
    kernel = c("gaussian", "matern52", "matern32"),
    loss = c("brier", "log_loss"),
    n_multi_start = NULL, theta = NULL
){
  # ---- Argument checking ----
  if (missing(X) || missing(Y)) {
    stop("Arguments 'X' and 'Y' must be provided.")
  }
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("'X' must be a numeric matrix or data frame.")
  }
  if (!is.numeric(as.matrix(X))) {
    stop("'X' must contain numeric values only.")
  }
  if (!is.matrix(Y) && !is.data.frame(Y)) {
    stop("'Y' must be a numeric matrix or data frame.")
  }
  if (!is.numeric(as.matrix(Y))) {
    stop("'Y' must contain numeric values only.")
  }

  X <- as.matrix(X)
  Y <- as.matrix(Y)

  d <- ncol(X)
  q <- ncol(Y)
  n <- nrow(X)

  if (nrow(Y) != n) {
    stop("Number of rows in 'Y' must match number of rows in 'X'.")
  }
  if (any(Y < 0)) stop("'Y' must be nonnegative counts or frequencies.")
  if (anyNA(X) || anyNA(Y)) stop("Missing values are not allowed in 'X' or 'Y'.")

  if (q < 2) {
    stop("'Y' must have at least two columns (multinomial outcomes).")
  }
  if (q == 2) {
    warning("For binary data, consider using the BKP model instead of DKP.")
  }

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

  if (!is.numeric(p0) || any(p0 < 0) || !all.equal(sum(p0), 1)) {
    stop("'p0' must be numeric, nonnegative, and sum to 1.")
  }

  if (prior == "fixed" && is.null(p0) && length(p0) != q) {
    stop("For DKP with prior = 'fixed', you must provide 'p0' with length equal to number of classes.")
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
      lower  = rep(-3, d), # relaxed lower bound
      upper  = rep(3, d),  # relaxed upper bound
      prior = prior, r0 = r0, p0 = p0,
      Xnorm = Xnorm, Y = Y,
      model = "DKP", loss = loss, kernel = kernel,
      control= list(trace=0))

    # ---- Extract optimal kernel parameters and loss ----
    best_index <- which.min(opt_res$value)
    gamma_opt  <- as.numeric(opt_res[best_index, 1:d])
    theta_opt  <- 10^gamma_opt
    loss_min   <- opt_res$value[best_index]
  }else{
    # ---- Use user-provided theta ----
    theta_opt <- theta
    loss_min <- loss_fun(gamma = log10(theta_opt), Xnorm = Xnorm, Y = Y,
                         prior = prior, r0 = r0, p0 = p0,
                         model = "DKP", loss = loss, kernel = kernel)
  }

  # ---- Compute kernel matrix at optimized hyperparameters ----
  K <- kernel_matrix(Xnorm, theta = theta_opt, kernel = kernel)

  # ---- Compute prior parameters (alpha0 and beta0) ----
  alpha0 <- get_prior(prior = prior, model = "DKP", r0 = r0, p0 = p0, Y = Y, K = K)

  # ---- Compute posterior parameters ----
  alpha_n <- alpha0 + as.matrix(K %*% Y)

  # ---- Construct and return the fitted model object ----
  DKP_model <- list(
    theta_opt = theta_opt, kernel = kernel,
    loss = loss, loss_min = loss_min,
    X = X, Xnorm = Xnorm, Xbounds = Xbounds, Y = Y,
    prior = prior, r0 = r0, p0 = p0,
    alpha0 = alpha0, alpha_n = alpha_n
  )
  class(DKP_model) <- "DKP"
  return(DKP_model)
}
