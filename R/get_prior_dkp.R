#' Construct Prior Parameters for Dirichlet Kernel Process (DKP)
#'
#' @description
#' Constructs prior Beta distribution parameters (alpha0 and beta0) at each location
#' based on user-specified prior type: noninformative, fixed, or data-adaptive.
#'
#' @param prior   Type of prior: one of `"noninformative"`, `"fixed"`, or `"adaptive"`.
#' @param r0      Global precision parameter (used in `"fixed"` and `"adaptive"`).
#' @param p0      Global prior mean used when prior = "fixed".
#' @param Y       A numeric vector of observed successes.
#' @param K       A precomputed kernel matrix (n x n), typically from `kernel_matrix()`.
#'
#' @return A list with a numeric vectors: `alpha0` (of length n).
#' @examples
#' get_prior_dkp(prior = "noninformative", Y = rep(1,3))
#' @export


get_prior_dkp <- function(prior = c("noninformative", "fixed", "adaptive"),
                      r0, p0, Y, K = NULL) {

  prior <- match.arg(prior)

  p <- ncol(Y)

  if (prior == "noninformative") {
    alpha0 <- as.vector(rep(1,p))
  } else if (prior == "fixed") {
    # Check input validity
    if (r0 <= 0) stop("r0 must be positive.")
    if (any(p0 < 0) || sum(p0) > 1) stop("sum(p0) must be in (0,1).")

    alpha0 <- as.vector(r0 * p0)
  } else if (prior == "adaptive") {
    # Check input validity
    if (r0 <= 0) stop("r0 must be positive.")
    if (is.null(Y) || is.null(K)) {
      stop("Y and K must be provided for adaptive prior.")
    }
    if (!is.matrix(K) || ncol(K) != nrow(Y)) stop("K must be a square matrix matching Y.")
    # Normalize K row-wise to get weights
    W <- K / rowSums(K)  # n x n matrix of weights
    p_hat <- W %*% (Y / rowSums(Y))   # estimated prior mean
    r_hat <- r0 * as.vector(rowSums(K)) + 1e-10   # estimated prior precision
    alpha0 <- r_hat * p_hat
  }
  return(list(alpha0 = alpha0))
}

