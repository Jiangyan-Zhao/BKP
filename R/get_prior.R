#' @title Construct Prior Parameters for BKP and DKP Models
#'
#' @description Construct prior shape parameters for Beta Kernel Process
#'   (BKP) and Dirichlet Kernel Process (DKP) models. The function supports
#'   \code{prior = "noninformative"}, \code{"fixed"}, and
#'   \code{"adaptive"} prior specifications.
#'
#' @inheritParams fit_BKP
#' @inheritParams fit_DKP
#' @param K Optional precomputed kernel matrix, typically obtained from
#'   \code{\link{kernel_matrix}}. For adaptive priors, \code{K} is required and
#'   must have one column per observed training input. The number of rows of
#'   \code{K} determines the number of locations at which prior parameters are
#'   evaluated. For noninformative and fixed priors, if \code{K = NULL}, the
#'   output length is inferred from \code{y}/\code{m} for BKP or from
#'   \code{Y} for DKP when available; otherwise a single prior row is returned.
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
#' - If \code{model = "DKP"}: a matrix \code{alpha0} of prior Dirichlet
#'   parameters at each input location (\code{n × q}).
#'
#' @details
#' The noninformative prior sets all Beta or Dirichlet shape parameters to 1.
#'
#' The fixed prior uses a common prior distribution at every evaluated
#' location. For BKP, this is a Beta prior with shape parameters
#' \code{r0 * p0} and \code{r0 * (1 - p0)}. For DKP, this is a Dirichlet
#' prior with concentration vector \code{r0 * p0}.
#'
#' The adaptive prior estimates location-specific prior means from the observed
#' data. BKP uses kernel smoothing of the observed proportions \code{y / m},
#' while DKP uses kernel smoothing of the row-normalized class counts in
#' \code{Y}. In both cases, the adaptive prior precision is scaled by the
#' local kernel mass.
#'
#' @references
#' Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#' Kernel Process Modeling}. arXiv.  <doi:10.48550/arXiv.2508.10447>.
#'
#' @examples
#' # -------------------------- BKP ---------------------------
#' set.seed(123)
#'
#' n <- 10
#' X <- matrix(runif(n * 2), ncol = 2)
#' y <- rbinom(n, size = 5, prob = 0.6)
#' m <- rep(5, n)
#' K <- kernel_matrix(X, theta = 0.2, kernel = "gaussian")
#'
#' # Adaptive BKP prior evaluated at the rows of K
#' prior_bkp <- get_prior(
#'   model = "BKP",
#'   prior = "adaptive",
#'   r0 = 2,
#'   y = y,
#'   m = m,
#'   K = K
#' )
#'
#' # Fixed BKP prior; with K = NULL, output length is inferred from y/m
#' prior_bkp_fixed <- get_prior(
#'   model = "BKP",
#'   prior = "fixed",
#'   r0 = 2,
#'   p0 = 0.5,
#'   y = y,
#'   m = m
#' )
#'
#' # -------------------------- DKP ---------------------------
#' n <- 15
#' q <- 3
#' X <- matrix(runif(n * 2), ncol = 2)
#'
#' true_pi <- t(apply(X, 1, function(x) {
#'   raw <- c(
#'     exp(-sum((x - 0.2)^2)),
#'     exp(-sum((x - 0.5)^2)),
#'     exp(-sum((x - 0.8)^2))
#'   )
#'   raw / sum(raw)
#' }))
#'
#' m <- sample(10:20, n, replace = TRUE)
#' Y <- t(sapply(seq_len(n), function(i) {
#'   rmultinom(1, size = m[i], prob = true_pi[i, ])
#' }))
#' K <- kernel_matrix(X, theta = 0.2, kernel = "gaussian")
#'
#' # Adaptive DKP prior evaluated at the rows of K
#' prior_dkp <- get_prior(
#'   model = "DKP",
#'   prior = "adaptive",
#'   r0 = 2,
#'   Y = Y,
#'   K = K
#' )
#'
#' # Fixed DKP prior; with K = NULL, output rows are inferred from Y
#' prior_dkp_fixed <- get_prior(
#'   model = "DKP",
#'   prior = "fixed",
#'   r0 = 2,
#'   p0 = rep(1 / q, q),
#'   Y = Y
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

  ## ---- Common argument checks ----
  if (!is.numeric(r0) || length(r0) != 1 || is.na(r0) || !is.finite(r0) || r0 <= 0) {
    stop("'r0' must be a positive finite scalar.")
  }

  if (!is.null(p0)) {
    if (!is.numeric(p0) || anyNA(p0) || any(!is.finite(p0)) || any(p0 < 0)) {
      stop("'p0' must be a numeric vector with nonnegative finite values.")
    }
  }

  if (!is.null(K)) {
    if (!is.matrix(K) || !is.numeric(K) || anyNA(K) || any(!is.finite(K))) {
      stop("'K' must be a numeric matrix with finite values and no NA values.")
    }
  }

  if (model == "BKP") {
    ## ---------------------- BKP checks ----------------------------------------
    if (!is.null(y)) {
      if (!is.numeric(y) || anyNA(y) || any(!is.finite(y))) {
        stop("'y' must be a numeric vector with finite values and no NA values.")
      }
      if (any(y < 0)) {
        stop("'y' must contain nonnegative counts.")
      }
    }
    if (!is.null(m)) {
      if (!is.numeric(m) || anyNA(m) || any(!is.finite(m))) {
        stop("'m' must be a numeric vector with finite values and no NA values.")
      }
      if (any(m <= 0)) {
        stop("'m' must contain positive trial counts.")
      }
    }
    if (!is.null(y) && !is.null(m)) {
      if (length(y) != length(m)) {
        stop("'y' and 'm' must have the same length.")
      }
      if (any(y > m)) {
        stop("'y' cannot exceed 'm'.")
      }
    }

    if (prior == "fixed") {
      if (is.null(p0) || !is.numeric(p0) || length(p0) != 1 ||
          is.na(p0) || !is.finite(p0) || p0 <= 0 || p0 >= 1) {
        stop("For fixed prior in BKP, 'p0' must be a scalar in (0, 1).")
      }
    }

    if (prior == "adaptive") {
      if (is.null(y) || is.null(m) || is.null(K)) {
        stop("For adaptive prior in BKP, 'y', 'm', and 'K' must be provided.")
      }
      if (ncol(K) != length(y)) {
        stop("'K' must have ncol = length(y).")
      }
    }

    if (is.null(K) && prior != "adaptive") {
      n_out <- 1L
      if (!is.null(y)) {
        n_out <- length(y)
      } else if (!is.null(m)) {
        n_out <- length(m)
      }
      K <- matrix(1, nrow = n_out, ncol = 1)
    }

    ## ----------------------- Call C++ computational core ----------------
    res <- get_prior_rcpp(
      model = model,
      prior = prior,
      r0 = r0,
      p0 = p0,
      y = y,
      m = m,
      Y = NULL,
      K = K
    )

    res$alpha0 <- as.vector(res$alpha0)
    res$beta0  <- as.vector(res$beta0)

    return(res)
  }else{
    ## ---------------------- DKP checks ----------------------------------------
    if (!is.null(Y)) {
      if (!is.matrix(Y) || !is.numeric(Y) || anyNA(Y) || any(!is.finite(Y))) {
        stop("'Y' must be a numeric matrix with finite values and no NA values.")
      }
      if (nrow(Y) < 1L || ncol(Y) < 2L) {
        stop("'Y' must have at least one row and at least two columns.")
      }
      if (any(Y < 0)) {
        stop("'Y' must contain nonnegative counts.")
      }
      if (any(rowSums(Y) <= 0)) {
        stop("Each row of 'Y' must have a positive row sum.")
      }
      q <- ncol(Y)
    } else if (!is.null(p0)) {
      q <- length(p0)
    } else {
      stop("Either 'Y' or 'p0' must be provided to determine the number of classes.")
    }

    if (q < 2L) {
      stop("For DKP, the number of classes must be at least two.")
    }

    if (prior == "fixed") {
      if (is.null(p0) || !is.numeric(p0) || length(p0) != q ||
          anyNA(p0) || any(!is.finite(p0)) ||
          any(p0 < 0) || abs(sum(p0) - 1) > 1e-10) {
        stop("For fixed prior in DKP, 'p0' must be a nonnegative finite numeric vector of length equal to the number of classes and sum to 1.")
      }
    }

    if (prior == "adaptive") {
      if (is.null(Y) || is.null(K)) {
        stop("'Y' and 'K' must be provided for adaptive prior in DKP.")
      }
      if (ncol(K) != nrow(Y)) {
        stop("'K' must have ncol = nrow(Y).")
      }
    }

    if (is.null(K) && prior != "adaptive") {
      n_out <- if (!is.null(Y)) nrow(Y) else 1L
      K <- matrix(1, nrow = n_out, ncol = 1)
    }

    ## ----------------------- Call C++ computational core ----------------
    res <- get_prior_rcpp(
      model = model,
      prior = prior,
      r0 = r0,
      p0 = p0,
      y = NULL,
      m = NULL,
      Y = Y,
      K = K
    )

    return(res$alpha0)
  }
}
