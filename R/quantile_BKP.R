#' @name quantile
#'
#' @title Posterior Quantiles from Fitted BKP Package Models
#'
#' @description Compute posterior quantiles from fitted \code{BKP},
#'   \code{DKP}, \code{TwinBKP}, and \code{TwinDKP} model objects. For
#'   \code{BKP} and \code{TwinBKP} objects, this returns posterior quantiles of
#'   the latent success probability. For \code{DKP} and \code{TwinDKP} objects,
#'   this returns marginal posterior quantiles for each class probability.
#'
#' @param x A fitted model object of class \code{"BKP"}, \code{"DKP"},
#'   \code{"TwinBKP"}, or \code{"TwinDKP"}, typically returned by
#'   \code{\link{fit_BKP}}, \code{\link{fit_DKP}},
#'   \code{\link{fit_TwinBKP}}, or \code{\link{fit_TwinDKP}}.
#' @param probs Numeric vector of probabilities specifying which posterior
#'   quantiles to return. Defaults to \code{c(0.025, 0.5, 0.975)}.
#' @param ... Additional arguments (currently unused).
#'
#' @return For \code{BKP} and \code{TwinBKP}: a numeric vector if
#'   \code{length(probs) = 1}, or a numeric matrix if
#'   \code{length(probs) > 1}, of posterior quantiles. Rows correspond to
#'   observations, and columns correspond to the requested probabilities.
#'
#'   For \code{DKP} and \code{TwinDKP}: a numeric matrix if
#'   \code{length(probs) = 1}, or a three-dimensional array if
#'   \code{length(probs) > 1}, of marginal posterior quantiles for class
#'   probabilities. Dimensions correspond to observations \eqn{\times} classes
#'   \eqn{\times} probabilities.
#'
#' @details For \code{BKP} and \code{TwinBKP} models, posterior quantiles are
#'   computed from the corresponding Beta posterior for the latent success
#'   probability. For \code{DKP} and \code{TwinDKP} models, marginal posterior
#'   quantiles for each class probability are computed from the Beta marginal
#'   distributions of the posterior Dirichlet distribution. These are marginal
#'   class-wise quantiles, not joint Dirichlet credible regions.
#'
#' @seealso \code{\link{fit_BKP}}, \code{\link{fit_DKP}},
#'   \code{\link{fit_TwinBKP}}, and \code{\link{fit_TwinDKP}} for model fitting;
#'   \code{\link{predict.BKP}}, \code{\link{predict.DKP}},
#'   \code{\link{predict.TwinBKP}}, and \code{\link{predict.TwinDKP}} for
#'   posterior prediction.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}. arXiv. <doi:10.48550/arXiv.2508.10447>.
#'
#' @keywords BKP
#'
#' @examples
#' # -------------------------- BKP and TwinBKP ---------------------------
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
#' # A fixed theta is used here only to keep the example fast and reproducible.
#' # In practice, omit theta to select it by leave-one-out cross-validation.
#' model <- fit_BKP(X, y, m, Xbounds = Xbounds, theta = 0.04)
#'
#' # Extract posterior quantiles
#' quantile(model)
#'
#' \dontrun{
#' # Larger TwinBKP example
#' n <- 1000
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit TwinBKP model using the default global lengthscale tuning
#' model <- fit_TwinBKP(
#'      X, y, m,
#'      Xbounds = Xbounds
#'    )
#'
#' # Extract posterior quantiles
#' quantile(model)
#' }
#'
#' @export
#' @method quantile BKP

quantile.BKP <- function(x, probs = c(0.025, 0.5, 0.975), ...) {
  # arguments checking
  if (!is.numeric(probs) || length(probs) < 1L ||
      anyNA(probs) || any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1)) {
    stop("'probs' must be a nonempty finite numeric vector with all values in [0, 1].")
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

