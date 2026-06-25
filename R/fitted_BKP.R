#' @name fitted
#'
#' @title Extract BKP, DKP, or TwinBKP Model Fitted Values
#'
#' @description Compute posterior fitted values from a fitted \code{BKP},
#'   \code{DKP}, or \code{TwinBKP} object. For \code{BKP} and
#'   \code{TwinBKP} objects, this returns the posterior mean probability of
#'   the positive class. For a \code{DKP} object, this returns the posterior
#'   mean probabilities for each class.
#'
#' @param object An object of class \code{BKP}, \code{DKP}, or
#'   \code{TwinBKP}, typically returned by \code{\link{fit_BKP}},
#'   \code{\link{fit_DKP}}, or \code{\link{fit_TwinBKP}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A numeric vector for \code{BKP} and \code{TwinBKP}, or a numeric
#'   matrix for \code{DKP}, containing posterior mean estimates at the
#'   training inputs.
#'
#' @details For \code{BKP} and \code{TwinBKP} models, the fitted values
#'   correspond to the posterior mean probability of the positive class,
#'   computed from the corresponding Beta posterior. For a \code{DKP} model,
#'   the fitted values correspond to the posterior mean probabilities for each
#'   class, derived from the posterior Dirichlet distribution of the class
#'   probabilities.
#'
#' @keywords BKP DKP TwinBKP
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv. \doi{10.48550/arXiv.2508.10447}
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
#' # Extract fitted values
#' fitted(model)
#'
#' @export
#' @method fitted BKP

fitted.BKP <- function(object, ...) {
  # Posterior beta parameters
  alpha_n <- object$alpha_n
  beta_n  <- object$beta_n

  # Posterior mean
  fitted_value <- alpha_n / (alpha_n + beta_n)

  return(fitted_value)
}

