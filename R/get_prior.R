#' @title Construct Prior Parameters for BKP/DKP Models
#'
#' @description Computes prior parameters for the Beta Kernel Process (BKP, for
#'   binary outcomes) or Dirichlet Kernel Process (DKP, for multi-class
#'   outcomes). Supports \code{prior = "noninformative"}, \code{"fixed"}, and
#'   \code{"adaptive"} strategies.
#'
#' @inheritParams fit_BKP
#' @inheritParams fit_DKP
#' @param K A precomputed kernel matrix, typically obtained from
#'   \code{\link{kernel_matrix}}. Can be rectangular (\code{m × n}), where
#'   \code{n} is the number of observed points and \code{m} the number of
#'   prediction locations.
#' @param model A character string specifying the model type: \code{"BKP"}
#'   (binary outcome) or \code{"DKP"} (multi-class outcome).
#' @param p0 For BKP, a scalar in \code{(0,1)} specifying the prior mean of
#'   success probability when \code{prior = "fixed"}. For DKP, a numeric vector
#'   of length equal to the number of classes specifying the global prior mean,
#'   which must sum to 1.
#' @param Y A numeric matrix of observed class counts (\code{n × q}), required
#'   only when \code{model = "DKP"}, where \code{n} is the number of
#'   observations and \code{q} the number of classes.
#'
#' @return
#' - If \code{model = "BKP"}: a list with
#'   \describe{
#'     \item{\code{alpha0}}{Vector of prior alpha parameters for the Beta
#'       distribution, length \code{n}.}
#'     \item{\code{beta0}}{Vector of prior beta parameters for the Beta
#'       distribution, length \code{n}.}
#'   }
#' - If \code{model = "DKP"}: a list containing
#'   \describe{
#'     \item{\code{alpha0}}{Matrix of prior Dirichlet parameters at each input
#'       location (\code{n × q}).}
#'   }
#'
#' @details
#' - \code{prior = "noninformative"}: flat prior; all parameters set to 1.
#' - \code{prior = "fixed"}:
#'   - BKP: uniform Beta prior \code{Beta(r0 * p0, r0 * (1 - p0))} across locations.
#'   - DKP: all rows of \code{alpha0} set to \code{r0 * p0}.
#' - \code{prior = "adaptive"}:
#'   - BKP: prior mean estimated at each location via kernel smoothing of observed
#' proportions \code{y/m}, with precision \code{r0}.
#'   - DKP: prior parameters computed by kernel-weighted smoothing of observed
#' class frequencies in \code{Y}, scaled by \code{r0}.
#'
#' @references
#' Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#' Kernel Process Modeling}. arXiv.
#' https://doi.org/10.48550/arXiv.2508.10447
#'
#' @examples
#' # -------------------------- BKP ---------------------------
#' set.seed(123)
#' n <- 10
#' X <- matrix(runif(n * 2), ncol = 2)
#' y <- rbinom(n, size = 5, prob = 0.6)
#' m <- rep(5, n)
#' K <- kernel_matrix(X)
#' prior_bkp <- get_prior(
#'   model = "BKP", prior = "adaptive", r0 = 2, y = y, m = m, K = K
#' )
#'
#' # -------------------------- DKP ---------------------------
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
#'   model = "DKP", prior = "adaptive", r0 = 2, Y = Y, K = K
#' )
#'
#' @seealso \code{\link{fit_BKP}} for fitting Beta Kernel Process models,
#'   \code{\link{fit_DKP}} for fitting Dirichlet Kernel Process models,
#'   \code{\link{predict.BKP}} and \code{\link{predict.DKP}} for making
#'   predictions, \code{\link{kernel_matrix}} for computing kernel matrices used
#'   in prior construction.
#'
#' @export
get_prior <- function(prior = c("noninformative", "fixed", "adaptive"),
                      model = c("BKP", "DKP"),
                      r0 = 2, p0 = NULL, y = NULL, m = NULL, Y = NULL, K = NULL)
{
  # ---- Argument checking ----
  model <- match.arg(model)
  prior <- match.arg(prior)

  if (!is.numeric(r0) || length(r0) != 1 || r0 <= 0) {
    stop("'r0' must be a positive scalar.")
  }

  if (!is.null(p0)) {
    if (!is.numeric(p0) || any(p0 < 0)) {
      stop("'p0' must be numeric and nonnegative.")
    }
  }

  if (model == "BKP") {
    if (!is.null(y) && (!is.numeric(y) || anyNA(y))) {
      stop("'y' must be a numeric vector with no NA values.")
    }
    if (!is.null(m) && (!is.numeric(m) || anyNA(m))) {
      stop("'m' must be a numeric vector with no NA values.")
    }
  }else{
    if (!is.null(Y) && (!is.matrix(Y) || anyNA(Y))) {
      stop("'Y' must be a numeric matrix with no NA values.")
    }
  }

  if (!is.null(K) && (!is.matrix(K) || anyNA(K))) {
    stop("'K' must be a numeric matrix with no NA values.")
  }

  if (model == "BKP") {
    # ============ Binary case ============
    if (prior == "noninformative") {
      # Assign uniform prior for each prediction location
      nrowK <- if (!is.null(K)) nrow(K) else 1
      alpha0 <- rep(1, nrowK)
      beta0  <- rep(1, nrowK)
    } else if (prior == "fixed") {
      if (!is.numeric(p0) || length(p0) != 1 || p0 <= 0 || p0 >= 1) {
        stop("For fixed prior in BKP, 'p0' must be in (0,1).")
      }
      nrowK <- if (!is.null(K)) nrow(K) else 1
      alpha0 <- rep(r0 * p0, nrowK)
      beta0  <- rep(r0 * (1 - p0), nrowK)
    } else if (prior == "adaptive") {
      if (is.null(y) || is.null(m) || is.null(K)) {
        stop("For adaptive prior in BKP, 'y', 'm', and 'K' must be provided.")
      }
      if (length(y) != length(m)) stop("'y' and 'm' must have the same length.")
      if (ncol(K) != length(y)) stop("'K' must have ncol = length(y).")

      # Row-normalized kernel weights
      W <- K / pmax(rowSums(K), 1e-6)    # m * n

      # Estimated mean and precision
      p_hat <- as.vector(W %*% (y / m))    # Estimated prior mean
      r_hat <- r0 * pmax(rowSums(K), 1e-3) # Estimated prior precision

      alpha0 <- r_hat * p_hat
      beta0  <- r_hat * (1 - p_hat)
      alpha0 <- pmax(alpha0, 1e-2) # Avoid numerical issues
      beta0 <- pmax(beta0, 1e-2)   # Avoid numerical issues
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
      stop("Either 'Y' or 'p0' must be provided to determine class dimension q.")
    }

    if (prior == "noninformative") {
      # Return a constant prior for all prediction points (say, m points)
      m <- if (!is.null(K)) nrow(K) else 1
      alpha0 <- matrix(1, nrow = m, ncol = q)
    } else if (prior == "fixed") {
      # Validate inputs
      if (is.null(p0)) stop("'p0' must be provided for fixed prior in DKP.")
      if (length(p0) != q) stop("length(p0) must match number of classes.")
      if (!isTRUE(all.equal(sum(p0), 1, tolerance = 1e-8))) {
        stop("'p0' must sum to 1.")
      }

      m <- if (!is.null(K)) nrow(K) else 1
      alpha0 <- matrix(rep(r0 * p0, each = m), nrow = m, byrow = TRUE)
    } else if (prior == "adaptive") {
      # Validate inputs
      if (is.null(Y) || is.null(K)) stop("'Y' and 'K' must be provided for adaptive prior in DKP.")
      if (ncol(K) != nrow(Y)) stop("'K' must have ncol = nrow(Y).")
      if (any(rowSums(Y) == 0)) stop("each row of Y must have positive sum for adaptive prior.")

      # Normalize kernel weights
      W <- K / pmax(rowSums(K), 1e-6)     # m * n

      # Estimate local class proportions
      Pi_hat <- W %*% (Y / rowSums(Y))    # m * q

      # Estimate local precision
      r_hat <- r0 * pmax(rowSums(K), 1e-3) # m * 1

      # Compute prior parameters
      alpha0 <- Pi_hat * r_hat
    }
    alpha0 <- pmax(alpha0, 1e-2)  # Avoid numerical issues
    return(alpha0)
  }
}
