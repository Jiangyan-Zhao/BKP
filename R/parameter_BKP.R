#' @name parameter
#'
#' @title Extract Model Parameters from a Fitted BKP, DKP, or TwinBKP Model
#'
#' @description Retrieve key model parameters from a fitted \code{BKP},
#'   \code{DKP}, or \code{TwinBKP} object. For a \code{BKP} model, this
#'   includes the optimized kernel hyperparameters and posterior Beta
#'   parameters. For a \code{DKP} model, this includes the kernel
#'   hyperparameters and posterior Dirichlet parameters. For a
#'   \code{TwinBKP} model, this additionally includes the global and local
#'   kernel parameters, selected global subset, and approximation controls.
#'
#' @param object An object of class \code{BKP}, \code{DKP}, or
#'   \code{TwinBKP}, typically returned by \code{\link{fit_BKP}},
#'   \code{\link{fit_DKP}}, or \code{\link{fit_TwinBKP}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A named list containing model parameters. Common entries include:
#' \itemize{
#'   \item \code{theta}: Estimated kernel hyperparameters.
#'   \item \code{alpha_n}: Posterior Dirichlet/Beta \eqn{\alpha} parameters.
#'   \item \code{beta_n}: Posterior Beta \eqn{\beta} parameters, returned for
#'     BKP and TwinBKP objects.
#' }
#' For \code{TwinBKP} objects, the returned list also includes
#' \code{theta_g}, \code{theta_l}, \code{global_kernel}, \code{local_kernel},
#' \code{global_indices}, and \code{control}.
#'
#' @keywords BKP DKP TwinBKP
#'
#' @seealso \code{\link{fit_BKP}}, \code{\link{fit_DKP}}, and
#'   \code{\link{fit_TwinBKP}} for model fitting.
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
