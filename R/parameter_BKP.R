#' @name parameter
#'
#' @title Extract Model Parameters from a Fitted BKP or DKP Model
#'
#' @description Retrieve the key model parameters from a fitted \code{BKP} or
#'   \code{DKP} object. For a \code{BKP} model, this typically includes the
#'   optimized kernel hyperparameters and posterior Beta parameters. For a
#'   \code{DKP} model, this includes the kernel hyperparameters and posterior
#'   Dirichlet parameters.
#'
#' @param object An object of class \code{BKP} or \code{DKP}, typically the
#'   result of a call to \code{\link{fit_BKP}} or \code{\link{fit_DKP}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{theta}: Estimated kernel hyperparameters.
#'   \item \code{alpha_n}: Posterior Dirichlet/Beta \eqn{\alpha} parameters.
#'   \item \code{beta_n}: (BKP only) Posterior Beta \eqn{\beta} parameters.
#' }
#'
#' @keywords BKP DKP
#'
#' @seealso \code{\link{fit_BKP}} for fitting BKP models, \code{\link{fit_DKP}}
#'   for fitting DKP models.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv.
#'   https://doi.org/10.48550/arXiv.2508.10447.
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
#' # Extract posterior and kernel parameters
#' parameter(model)
#'
#' @export
parameter <- function(object, ...) {
  UseMethod("parameter")
}

#' @rdname parameter
#' @export
#' @method parameter BKP
parameter.BKP <- function(object, ...) {
  list(
    theta   = object$theta_opt,
    alpha_n = object$alpha_n,
    beta_n  = object$beta_n
  )
}
