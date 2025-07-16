#' @title Construct Prior Parameters for Dirichlet Kernel Process (DKP)
#'
#' @description
#' Constructs prior Dirichlet distribution parameters (alpha0) at each location
#' based on user-specified prior type: noninformative, fixed, or data-adaptive.
#'
#' @param prior   Type of prior: one of `"noninformative"`, `"fixed"`, or `"adaptive"`.
#' @param r0      Global precision parameter (used in `"fixed"` and `"adaptive"`).
#' @param p0      Global prior mean used when prior = "fixed".
#' @param Y       A numeric vector of observed successes.
#' @param K       A precomputed kernel matrix (n x n), typically from `kernel_matrix()`.
#'
#' @return A list with a numeric vectors: `alpha0` (of length n).
#'
#' @examples
#' # Simulated multi-class data
#' set.seed(123)
#' n <- 15           # number of training points
#' p <- 3            # number of classes
#' X <- matrix(runif(n * 2), ncol = 2)
#'
#' # Simulate class probabilities and draw multinomial counts
#' true_pi <- t(apply(X, 1, function(x) {
#'   raw <- c(
#'     exp(-sum((x - 0.2)^2)),
#'     exp(-sum((x - 0.5)^2)),
#'     exp(-sum((x - 0.8)^2))
#'   )
#'   raw / sum(raw)
#' }))
#' m <- sample(10:20, n, replace = TRUE)
#' Y <- t(apply(true_pi, 1, function(p) rmultinom(1, size = sample(10:20, 1), prob = p)))
#'
#' # Compute kernel matrix (e.g., Gaussian)
#' K <- kernel_matrix(X, theta = rep(0.2, 2), kernel = "gaussian")
#'
#' # Construct prior (adaptive)
#' prior_dkp <- get_prior_dkp(prior = "adaptive", r0 = 2, Y = Y, K = K)
#' str(prior_dkp)
#'
#' @export


get_prior_dkp <- function(prior = c("noninformative", "fixed", "adaptive"),
                      r0 = 2, p0 = NULL, Y = NULL, K = NULL) {

  prior <- match.arg(prior)

  if (!is.null(Y)) {
    Y <- as.matrix(Y)
    q <- ncol(Y)
  } else if (!is.null(p0)) {
    q <- length(p0)
  } else {
    stop("Either Y or p0 must be provided to determine class dimension q.")
  }


  if (prior == "noninformative") {
    # Return a constant prior for all prediction points (say, m points)
    m <- if (!is.null(K)) nrow(K) else 1
    alpha0 <- matrix(1, nrow = m, ncol = q)

  } else if (prior == "fixed") {
    # Validate inputs
    if (r0 <= 0) stop("r0 must be positive.")
    if (is.null(p0)) stop("p0 must be provided for fixed prior.")
    if (length(p0) != q) stop("Length of p0 must match the number of classes.")
    if (any(p0 < 0) || abs(sum(p0) - 1) > 1e-6) {
      stop("p0 must be a valid probability vector (nonnegative, sums to 1).")
    }

    m <- if (!is.null(K)) nrow(K) else 1
    alpha0 <- matrix(rep(r0 * p0, each = m), nrow = m, byrow = TRUE)
  } else if (prior == "adaptive") {
    # Validate inputs
    if (is.null(Y) || is.null(K)) stop("Y and K must be provided for adaptive prior.")
    if (!is.matrix(K) || ncol(K) != nrow(Y)) {
      stop("K must be an m * n matrix where n = nrow(Y).")
    }
    if (r0 <= 0) stop("r0 must be positive.")

    # Normalize kernel weights
    W <- K / rowSums(K)  # m * n

    # Estimate local class proportions
    Pi_hat <- W %*% (Y / rowSums(Y))  # m * q

    # Estimate local precision
    r_hat <- r0 * rowSums(K) + 1e-10  # m * 1

    # Compute prior parameters
    alpha0 <- Pi_hat * r_hat
  }
  return(alpha0 = alpha0)
}

