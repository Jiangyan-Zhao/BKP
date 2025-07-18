#' @title Construct Prior Parameters for the BKP Model
#'
#' @description Computes the prior Beta distribution parameters \code{alpha0}
#'   and \code{beta0} at each input location, based on the chosen prior
#'   specification. Supports noninformative, fixed, and data-adaptive prior
#'   strategies.
#'
#' @param prior Character string specifying the type of prior to use. One of
#'   \code{"noninformative"}, \code{"fixed"}, or \code{"adaptive"}.
#' @param r0 Positive scalar indicating the global precision parameter. Used
#'   when \code{prior} is \code{"fixed"} or \code{"adaptive"}.
#' @param p0 Prior mean for the success probability (in (0,1)). Used only when
#'   \code{prior = "fixed"}.
#' @param y Numeric vector of observed successes, of length \code{n}.
#' @param m Numeric vector of total binomial trials, of length \code{n}.
#' @param K A precomputed kernel matrix of size \code{n × n}, typically obtained
#'   from \code{\link{kernel_matrix}}.
#'
#' @return A list with two numeric vectors:
#' \describe{
#'   \item{\code{alpha0}}{Prior alpha parameters of the Beta distribution, length \code{n}.}
#'   \item{\code{beta0}}{Prior beta parameters of the Beta distribution, length \code{n}.}
#' }
#'
#' @details
#' - For \code{prior = "noninformative"}, all prior parameters are set to 1 (noninformative prior).
#' - For \code{prior = "fixed"}, all locations share the same Beta prior:
#' \code{Beta(r0 * p0, r0 * (1 - p0))}.
#' - For \code{prior = "adaptive"}, the prior mean at each location is computed by
#' kernel smoothing the observed proportions \code{y/m}, and precision \code{r0}
#' is distributed accordingly.
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
#' # Construct adaptive prior
#' prior <- get_prior(prior = "adaptive", r0 = 2, y = y, m = m, K = K)
#'
#' @seealso \code{\link{get_prior_dkp}}, \code{\link{fit.BKP}},
#'   \code{\link{predict.BKP}}, \code{\link{kernel_matrix}}
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
