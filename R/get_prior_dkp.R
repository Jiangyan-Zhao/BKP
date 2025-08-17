#' @title Construct Prior Parameters for the DKP Model
#'
#' @description Computes prior Dirichlet distribution parameters \code{alpha0}
#'   at each input location for the Dirichlet Kernel Process model, based on the
#'   specified prior type: noninformative, fixed, or adaptive.
#'
#' @inheritParams get_prior
#' @param p0 Numeric vector specifying the global prior mean for each class
#'   (must sum to 1). Only used when \code{prior = "fixed"}. Should be of length
#'   equal to the number of classes.
#' @param Y Numeric matrix of observed class counts of size \code{n × q}, where
#'   \code{n} is the number of observations and \code{q} the number of classes.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{alpha0}}{A numeric matrix of prior Dirichlet parameters at each input location;
#'     dimension \code{n × q}.}
#' }
#'
#' @details
#' - When \code{prior = "noninformative"}, all entries in \code{alpha0} are set to 1 (flat Dirichlet).
#' - When \code{prior = "fixed"}, all rows of \code{alpha0} are set to \code{r0 * p0}.
#' - When \code{prior = "adaptive"}, each row of \code{alpha0} is computed by kernel-weighted
#' smoothing of the observed relative frequencies in \code{Y}, scaled by
#' \code{r0}.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv.
#'   https://doi.org/10.48550/arXiv.2508.10447.
#'
#' @examples
#' # Simulated multi-class data
#' set.seed(123)
#' n <- 15           # number of training points
#' q <- 3            # number of classes
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
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Compute kernel matrix (Gaussian)
#' K <- kernel_matrix(X, theta = rep(0.2, 2), kernel = "gaussian")
#'
#' # Construct adaptive prior
#' prior_dkp <- get_prior_dkp(prior = "adaptive", r0 = 2, Y = Y, K = K)
#'
#' @seealso \code{\link{get_prior}}, \code{\link{fit.DKP}},
#'   \code{\link{predict.DKP}}, \code{\link{kernel_matrix}}
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
    r_hat <- r0 * rowSums(K)          # m * 1

    # Compute prior parameters
    alpha0 <- Pi_hat * r_hat
  }
  alpha0 <- pmax(alpha0, 1e-3)
  return(alpha0 = alpha0)
}

