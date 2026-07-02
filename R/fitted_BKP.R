#' @name fitted
#'
#' @title Extract Fitted Posterior Means from BKP Package Model Objects
#'
#' @description Extract posterior fitted values from fitted \code{BKP},
#'   \code{DKP}, \code{TwinBKP}, or \code{TwinDKP} model objects. For
#'   \code{BKP} and \code{TwinBKP} objects, this returns the posterior mean
#'   success probability at each training input. For \code{DKP} and
#'   \code{TwinDKP} objects, this returns the posterior mean class-probability
#'   vector at each training input.
#'
#' @param object A fitted model object of class \code{BKP}, \code{DKP},
#'   \code{TwinBKP}, or \code{TwinDKP}, typically returned by
#'   \code{\link{fit_BKP}}, \code{\link{fit_DKP}},
#'   \code{\link{fit_TwinBKP}}, or \code{\link{fit_TwinDKP}}.
#' @param ... Additional arguments. Currently unused.
#'
#' @return For \code{BKP} and \code{TwinBKP} objects, a numeric vector
#'   containing posterior mean success probabilities at the training inputs.
#'   For \code{DKP} and \code{TwinDKP} objects, a numeric matrix containing
#'   posterior mean class probabilities at the training inputs, with one row per
#'   training input and one column per class.
#'
#' @details The \code{fitted()} method extracts posterior means already stored
#'   in the fitted model object. For \code{BKP} and \code{TwinBKP}, the fitted
#'   value at a training input is computed from the corresponding Beta posterior
#'   shape parameters. For \code{DKP} and \code{TwinDKP}, fitted values are
#'   obtained by normalizing the corresponding Dirichlet posterior concentration
#'   parameters across classes.
#'
#' @seealso \code{\link{fit_BKP}}, \code{\link{fit_DKP}},
#'   \code{\link{fit_TwinBKP}}, and \code{\link{fit_TwinDKP}} for model
#'   fitting; \code{\link{predict}} for posterior prediction at new input locations.
#'
#' @keywords BKP
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv. <doi:10.48550/arXiv.2508.10447>.
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
#' # Extract fitted values
#' fitted(model)
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
#' # Extract fitted values
#' fitted(model)
#' }
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

