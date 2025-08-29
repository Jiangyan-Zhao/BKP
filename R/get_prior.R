#' @title Construct Prior Parameters for BKP/DKP Models
#'
#' @description Computes the prior parameters for the Beta Kernel Process (BKP,
#'   binary outcomes) or Dirichlet Kernel Process (DKP, multi-class outcomes).
#'   The function supports \code{prior = "noninformative"}, \code{"fixed"}, and
#'   \code{"adaptive"} strategies.
#'
#' @inheritParams fit.BKP
#' @inheritParams fit.DKP
#' @param K A precomputed kernel matrix of size \code{n × n}, typically obtained
#'   from \code{\link{kernel_matrix}}.
#' @param model_type A character string, either \code{"BKP"} (binary outcome) or
#'   \code{"DKP"} (multi-class outcome).
#' @param p0
#'   - For BKP: a scalar in \code{(0,1)} specifying the prior mean of success
#'   probability when \code{prior = "fixed"}.
#'   - For DKP: a numeric vector of length equal to the number of classes,
#'   specifying the global prior mean (must sum to 1).
#' @param Y A numeric matrix of observed class counts of size \code{n × q} (only
#'   required when \code{model_type = "DKP"}), where \code{n} is the number of
#'   observations and \code{q} the number of classes.
#'
#' @return
#' - If \code{model_type = "BKP"}: a list with two numeric vectors
#'   \describe{
#'     \item{\code{alpha0}}{Prior alpha parameters of the Beta distribution,
#'       length \code{n}.}
#'     \item{\code{beta0}}{Prior beta parameters of the Beta distribution,
#'       length \code{n}.}
#'   }
#' - If \code{model_type = "DKP"}: a list containing
#'   \describe{
#'     \item{\code{alpha0}}{A numeric matrix of prior Dirichlet parameters at
#'       each input location, dimension \code{n × q}.}
#'   }
#'
#' @details
#' - When \code{prior = "noninformative"}:
#'   - BKP: all prior parameters are set to 1 (flat Beta).
#'   - DKP: all entries in \code{alpha0} are set to 1 (flat Dirichlet).
#' - When \code{prior = "fixed"}:
#'   - BKP: all locations share the same Beta prior
#' \code{Beta(r0 * p0, r0 * (1 - p0))}.
#'   - DKP: all rows of \code{alpha0} are set to \code{r0 * p0}.
#' - When \code{prior = "adaptive"}:
#'   - BKP: the prior mean at each location is computed by kernel smoothing
#' of the observed proportions \code{y/m}, with precision \code{r0}.
#'   - DKP: each row of \code{alpha0} is computed by kernel-weighted smoothing
#' of the observed relative frequencies in \code{Y}, scaled by \code{r0}.
#'
#' @references
#' Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#' Kernel Process Modeling}. arXiv.
#' https://doi.org/10.48550/arXiv.2508.10447
#'
#' @examples
#' ## BKP example
#' set.seed(123)
#' n <- 10
#' X <- matrix(runif(n * 2), ncol = 2)
#' y <- rbinom(n, size = 5, prob = 0.6)
#' m <- rep(5, n)
#' K <- kernel_matrix(X)
#' prior_bkp <- get_prior(
#'   model_type = "BKP", prior = "adaptive", r0 = 2, y = y, m = m, K = K
#' )
#'
#' ## DKP example
#' set.seed(123)
#' n <- 15; q <- 3
#' X <- matrix(runif(n * 2), ncol = 2)
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
#' K <- kernel_matrix(X, theta = rep(0.2, 2), kernel = "gaussian")
#' prior_dkp <- get_prior(
#'   model_type = "DKP", prior = "adaptive", r0 = 2, Y = Y, K = K
#' )
#'
#' @seealso \code{\link{fit.BKP}}, \code{\link{fit.DKP}},
#'   \code{\link{predict.BKP}}, \code{\link{predict.DKP}},
#'   \code{\link{kernel_matrix}}
#'
#' @export
get_prior <- function(prior = c("noninformative", "fixed", "adaptive"),
                      model_type = c("BKP", "DKP"),
                      r0 = 2, p0 = 0.5, y = NULL, m = NULL, Y = NULL, K = NULL)
{
  model_type <- match.arg(model_type)
  prior <- match.arg(prior)

  if (model_type == "BKP") {
    # ============ Binary case ============
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
      p_hat <- as.vector(W %*% (y / m))   # Estimated prior mean
      r_hat <- r0 * rowSums(K)            # Estimated prior precision

      alpha0 <- r_hat * p_hat
      beta0  <- r_hat * (1 - p_hat)
      alpha0 <- pmax(alpha0, 1e-3)
      beta0 <- pmax(beta0, 1e-3)
    }
    return(list(alpha0 = alpha0, beta0 = beta0))
  } else {
    # ============ Multiclass case ============
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
}
