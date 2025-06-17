#' @title Fit a Beta Kernel Process (BKP) model
#'
#' @description
#' This function fits a Beta Kernel Process (BKP) model by optimizing kernel parameters
#' to minimize the Brier score. It supports different kernel types and allows for
#' flexible input of covariates (X), successes (y), and trials (m).
#'
#' @param data A data frame or matrix containing X, y, and m in order.
#'   If provided, individual `X`, `y`, and `m` parameters will be ignored.
#' @param X A data frame, matrix, or vector of covariates, with columns representing
#'   the dimensions of X. Only used if `data` is not provided.
#' @param y A sequence of numbers of successes out of `m` trials.
#'   Only used if `data` is not provided.
#' @param m A sequence of numbers of trials. Only used if `data` is not provided.
#' @param kernel.type Character string specifying the type of kernel to use.
#'   Accepted values are "gaussian", "matern52", or "matern32". Defaults to "gaussian".
#' @param loss.type Character string specifying the loss fucntion to use.
#'   Accepted values are "likelihood" or "brier". Defaults to "brier".
#' @param alpha0 A scalar or a sequence, the prior parameter for the beta distribution (shape1).
#'   Default by 1. If a sequence, its length must match the number of observations.
#' @param beta0 A scalar or a sequence, the prior parameter for the beta distribution (shape2).
#'   Default by 1. If a sequence, its length must match the number of observations.
#'
#' @return A `bkp` object, which is a list containing the optimized model parameters and
#'   input data:
#' \itemize{
#'   \item \code{bestTheta}: The optimized kernel parameters.
#'   \item \code{minLoss}: The minimum Brier score achieved during optimization.
#'   \item \code{kfun}: The selected kernel function based on `kernel.type`.
#'   \item \code{kernel.type}: The type of kernel used ("gaussian", "matern52", or "matern32").
#'   \item \code{y}: The input `y` values (number of successes).
#'   \item \code{m}: The input `m` values (number of trials).
#'   \item \code{normalizedX}: The covariate matrix `X` after min-max normalization to [0, 1].
#'   \item \code{originalX}: The original, unnormalized covariate matrix `X`.
#'   \item \code{alpha0}: The `alpha0` prior parameter(s) used for the beta distribution.
#'   \item \code{beta0}: The `beta0` prior parameter(s) used for the beta distribution.
#' }
#' @export
#' @importFrom optimx multistart
#' @importFrom lhs randomLHS
#'
#' @examples
#' ### 1D
#' set.seed(123)
#' n <- 100
#' x <- seq(-2, 2, length = n)
#' true_pi <- (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#' df <- data.frame(x = x, y = y, m = m)
#' xx <- seq(-2, 2, length = 100) #new data points
#' model <- fit.bkp(df)
#' head(predict(model,xx))
#' plot(model)
#' print(model)
#' summary(model)
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
#' x <- lhs::randomLHS(n = n, k = 2)
#' true_pi <- pnorm(f(x))
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#' df <- data.frame(x = x, y = y, m = m)
#' xx1 <- seq(0, 1, length.out = 100)
#' xx2 <- seq(0, 1, length.out = 100)
#' xx <- expand.grid(xx1 = xx1, xx2 = xx2)
#' model <- fit.bkp(df)
#' head(predict(model,xx))
#' plot(model)
#' print(model)
#' summary(model)
#'
#' #plot the true probability
#' true_pi <- pnorm(f(xx))
#' pred_prob_matrix <- matrix(true_pi,
#'                            nrow = length(xx1),
#'                            ncol = length(xx2))
#' graphics::image(xx1, xx2, pred_prob_matrix,
#'                 xlab ="X1", # Use X column name if available.
#'                 ylab ="X2", # Use X column name if available.
#'                 main = "True Probability",
#'                 col = grDevices::hcl.colors(100, "viridis"))

fit.bkp <- function(
    data,
    X = NULL,
    y = NULL,
    m = NULL,
    kernel.type = "gaussian",
    loss.type = "brier",
    alpha0 = 1,
    beta0 = 1
){
  # Handle input data: prioritize 'data' data frame, otherwise use individual X, y, m.
  # This block ensures that X, y, and m are correctly extracted regardless of input method.
  if (!missing(data)) { # Check if 'data' argument was explicitly provided by the user.
    nCol <- ncol(data)
    # Ensure 'data' has at least 3 columns (for X, y, and m).
    if (nCol < 3) stop("The 'data' frame must contain at least three columns (covariates X, y, and m in order).")
    d <- nCol - 2 # Determine the number of covariate dimensions (X columns).
    # Extract y (last column), m (second to last column), and X (remaining columns).
    y <- as.matrix(data[, d + 1])
    m <- as.matrix(data[, d + 2])
    X <- as.matrix(data[, 1:d])
  } else { # If 'data' is not provided, ensure X, y, and m are all supplied individually.
    if (is.null(X) || is.null(y) || is.null(m)) {
      stop("Either the 'data' argument must be supplied, or 'X', 'y', and 'm' must all be provided separately.")
    }
    # Convert individual X, y, m inputs to matrices for consistent processing.
    X <- as.matrix(X)
    d <- ncol(X) # Determine the number of covariate dimensions from X.
    y <- as.matrix(y)
    m <- as.matrix(m)
  }

  # Perform dimension consistency checks for X, y, and m.
  # All inputs must have the same number of observations.
  if (nrow(X) != length(y) || nrow(X) != length(m) || length(y) != length(m)) {
    stop("The sample sizes of X, y, and m do not match. Please ensure they have the same number of observations.")
  }
  # Check consistency for alpha0 and beta0: they must be scalars or vectors matching sample size.
  if ((length(alpha0) != 1 || length(beta0) != 1) & (length(alpha0) != length(m) || length(beta0) != length(m))) {
    stop("alpha0 and beta0 must be scalars or sequences with the same sample size as X, y, and m.")
  }

  n <- nrow(X) # Number of observations.
  originalX <- X # Store the original X before normalization.
  X <- normalize01(X) # Normalize X covariates to the [0, 1] range.
  # (Assuming normalize01 is an internal helper function)

  # Validate and select the kernel type.
  # match.arg ensures only accepted values are used and provides partial matching.
  kernel.type <- match.arg(kernel.type, c("gaussian", "matern52", "matern32"))

  # Define the specific kernel function based on 'kernel.type'.
  # (Assuming gaussianKernel, matern52Kernel, matern32Kernel are internal helper functions)
  kernelFunc <- switch(kernel.type,
                       gaussian = gaussianKernel,
                       matern52 = matern52Kernel,
                       matern32 = matern32Kernel
  )
  # Create a wrapper 'kfun' for the selected kernel, handling distance calculation.
  # (Assuming distFun is an internal helper function for distance calculation)
  kfun <- function(x1, x2, theta) { kernelFunc(distFun(x1, x2, theta)) }

  # Generate initial guesses for gamma parameters using Latin Hypercube Sampling.
  # These 'gamma' values are related to kernel parameters (theta = 10^(-gamma)).
  # The range of initialGamma is chosen to provide a good search space for optimization.
  initialGamma <- lhs::randomLHS(10 * d, d) * (2 + log10(500)) - 2 - log10(d)

  # Perform multi-start optimization to find the best kernel parameters.
  # 'optimx::multistart' explores multiple starting points to find a global optimum.
  # 'fn' is the objective function (Brier score) to be minimized.
  # 'brierScore' is an assumed internal function calculating the Brier score.
  # 'L-BFGS-B' is an optimization method that allows for box constraints (lower/upper bounds).
  if (loss.type == "likelihood"){
    res <- optimx::multistart(parmat = initialGamma,
                              fn = function(gamma) { logPseudoLikelihood(gamma, alpha0, beta0, kfun, X, y, m, n) },
                              method = "L-BFGS-B",
                              lower = -2 - log10(d), # Lower bound for gamma parameters.
                              upper = log10(500) - log10(d)) # Upper bound for gamma parameters.
  } else if (loss.type == "brier"){
    res <- optimx::multistart(parmat = initialGamma,
                              fn = function(gamma) { brierScore(gamma, alpha0, beta0, kfun, X, y, m, n) },
                              method = "L-BFGS-B",
                              lower = -2 - log10(d), # Lower bound for gamma parameters.
                              upper = log10(500) - log10(d)) # Upper bound for gamma parameters.
  }

  # Extract the results from the optimization.
  bestIndex <- which.min(res$value) # Find the index of the minimum loss.
  bestGamma <- as.numeric(res[bestIndex, 1:d]) # Get the gamma parameters corresponding to min loss.
  bestTheta <- 10^(-bestGamma) # Transform gamma back to the kernel parameters (theta).
  minLoss <- res$value[bestIndex] # Get the minimum loss value.

  # Construct the 'bkp' model object as a list.
  # This list contains all essential information about the fitted model.
  bkp.model <- list(
    bestTheta = bestTheta,
    minLoss = minLoss,
    kfun = kfun,
    kernel.type = kernel.type,
    loss.type = loss.type,
    y = y,
    m = m,
    normalizedX = X, # Return the normalized covariates.
    originalX = originalX, # Return the original covariates.
    alpha0 = alpha0,
    beta0 = beta0
  )
  # Assign the "bkp" class to the object, enabling S3 generic methods.
  class(bkp.model) <- "bkp"
  return(bkp.model) # Return the fitted BKP model object.
}
