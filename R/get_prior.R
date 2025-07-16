#' @title Construct Prior Parameters for Beta Kernel Process (BKP)
#'
#' @description Constructs prior Beta distribution parameters (alpha0 and beta0)
#' at each location based on user-specified prior type: noninformative, fixed,
#' or data-adaptive.
#'
#' @param prior   Type of prior: one of `"noninformative"`, `"fixed"`, or
#'   `"adaptive"`.
#' @param r0      Global precision parameter (used in `"fixed"` and
#'   `"adaptive"`).
#' @param p0      Global prior mean (used in `"fixed"` prior only; must be in
#'   (0,1)).
#' @param y       A numeric vector of observed successes.
#' @param m       A numeric vector of total trials at each location.
#' @param K       A precomputed kernel matrix (m x n), typically from
#'   `kernel_matrix()`.
#'
#' @return A list with two numeric vectors: `alpha0` and `beta0` (each of length
#'   m).
#'
#' @examples
#' # Simulated data
#' set.seed(123)
#' n <- 10
#' X <- matrix(runif(n * 2), ncol = 2)
#' y <- rbinom(n, size = 5, prob = 0.6)
#' m <- rep(5, n)
#'
#' # Example kernel matrix (Gaussian)
#' K <- kernel_matrix(X)
#'
#' # Construct prior (adaptive)
#' prior <- get_prior(prior = "adaptive", r0 = 2, y = y, m = m, K = K)
#' str(prior)
#'
#' @export

get_prior <- function(prior = c("noninformative", "fixed", "adaptive"),
                      r0 = 2, p0 = 0.5, y = NULL, m = NULL, K = NULL) {

  prior <- match.arg(prior)

  if (prior == "noninformative") {
    # Assign uniform prior for each prediction location
    nrowK <- if (!is.null(K)) nrow(K) else 1
    alpha0 <- rep(1, nrowK)
    beta0  <- rep(1, nrowK)
  } else if (prior == "fixed") {
    # Validate inputs
    if (r0 <= 0) stop("r0 must be positive.")
    if (p0 <= 0 || p0 >= 1) stop("p0 must be in (0, 1).")

    nrowK <- if (!is.null(K)) nrow(K) else 1
    alpha0 <- rep(r0 * p0, nrowK)
    beta0  <- rep(r0 * (1 - p0), nrowK)
  } else if (prior == "adaptive") {
    # Validate required inputs
    if (is.null(y) || is.null(m) || is.null(K)) {
      stop("y, m, and K must be provided for adaptive prior.")
    }
    if (length(y) != length(m)) stop("y and m must have the same length.")
    if (!is.matrix(K) || ncol(K) != length(y)) {
      stop("K must be an m * n matrix with n = length(y).")
    }
    if (r0 <= 0) stop("r0 must be positive.")

    # Row-normalized kernel weights
    W <- K / rowSums(K)   # m * n

    # Estimated mean and precision
    p_hat <- as.vector(W %*% (y / m))          # Estimated prior mean
    r_hat <- r0 * rowSums(K) + 1e-10           # Estimated prior precision

    alpha0 <- r_hat * p_hat
    beta0  <- r_hat * (1 - p_hat)
  }
  return(list(alpha0 = alpha0, beta0 = beta0))
}
