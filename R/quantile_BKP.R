#' @name quantile
#'
#' @title Posterior Quantiles from a Fitted BKP, DKP, or TwinBKP Model
#'
#' @description Compute posterior quantiles from a fitted \code{BKP},
#'   \code{DKP}, or \code{TwinBKP} model. For \code{BKP} and
#'   \code{TwinBKP} objects, this returns posterior quantiles of the positive
#'   class probability. For a \code{DKP} object, this returns posterior
#'   quantiles for each class probability.
#'
#' @param x An object of class \code{BKP}, \code{DKP}, or \code{TwinBKP},
#'   typically returned by \code{\link{fit_BKP}}, \code{\link{fit_DKP}}, or
#'   \code{\link{fit_TwinBKP}}.
#' @param probs Numeric vector of probabilities specifying which posterior
#'   quantiles to return. Defaults to \code{c(0.025, 0.5, 0.975)}.
#' @param ... Additional arguments (currently unused).
#'
#' @return For \code{BKP} and \code{TwinBKP}: a numeric vector if
#'   \code{length(probs) = 1}, or a numeric matrix if
#'   \code{length(probs) > 1}, of posterior quantiles. Rows correspond to
#'   observations, and columns correspond to the requested probabilities.
#'
#'   For \code{DKP}: a numeric matrix if \code{length(probs) = 1}, or a 3D
#'   array if \code{length(probs) > 1}, of posterior quantiles. Dimensions
#'   correspond to observations \eqn{\times} classes \eqn{\times}
#'   probabilities.
#'
#' @details For \code{BKP} and \code{TwinBKP} models, posterior quantiles are
#'   computed from the corresponding Beta posterior for the positive class
#'   probability. For a \code{DKP} model, marginal posterior quantiles for each
#'   class are computed from the exact Beta marginal distributions of the
#'   posterior Dirichlet distribution.
#'
#' @seealso \code{\link{fit_BKP}}, \code{\link{fit_DKP}}, and
#'   \code{\link{fit_TwinBKP}} for model fitting.
#'
#' @keywords BKP DKP TwinBKP
#'
#' @examples
#' # -------------------------- BKP ---------------------------
#' set.seed(123)
#'
#' # Define true success probability function
#' true_pi_fun <- function(x) {
#'   (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' }
#'
#' n <- 30
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit BKP model
#' model <- fit_BKP(X, y, m, Xbounds = Xbounds)
#'
#' # Extract posterior quantiles
#' quantile(model)
#'
#' @export
#' @method quantile BKP

quantile.BKP <- function(x, probs = c(0.025, 0.5, 0.975), ...) {
  # arguments checking
  if (!is.numeric(probs) || anyNA(probs) || any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1)) {
    stop("'probs' must be a finite numeric vector with all values in [0, 1].")
  }

  # Extract posterior beta parameters
  alpha_n <- x$alpha_n
  beta_n  <- x$beta_n

  n <- length(alpha_n)

  if (length(probs) > 1) {
    # Posterior quantiles matrix: rows = observations, cols = probs.
    post_q <- matrix(
      qbeta(
        rep(probs, each = n),
        rep(alpha_n, times = length(probs)),
        rep(beta_n, times = length(probs))
      ),
      nrow = n,
      ncol = length(probs)
    )
    colnames(post_q) <- paste0(probs * 100, "%")
  } else {
    # Single probability: return a vector
    post_q <- qbeta(probs, alpha_n, beta_n)
  }

  return(post_q)
}

