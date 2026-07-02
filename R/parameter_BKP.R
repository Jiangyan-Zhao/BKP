#' @name parameter
#'
#' @title Extract Model Parameters from Fitted BKP Package Models
#'
#' @description Extract key fitted parameters from \code{BKP}, \code{DKP},
#'   \code{TwinBKP}, and \code{TwinDKP} model objects. For full BKP and DKP
#'   models, this includes the fitted kernel lengthscale parameter(s) and
#'   posterior shape or concentration parameters. For TwinBKP and TwinDKP models,
#'   the returned list also includes global-local approximation parameters,
#'   selected global subset indices, and fitting controls.
#'
#' @param object A fitted model object of class \code{BKP}, \code{DKP},
#'   \code{TwinBKP}, or \code{TwinDKP}, typically returned by
#'   \code{\link{fit_BKP}}, \code{\link{fit_DKP}},
#'   \code{\link{fit_TwinBKP}}, or \code{\link{fit_TwinDKP}}.
#' @param ... Additional arguments. Currently unused.
#'
#' @return A named list containing fitted model parameters. Common entries are:
#' \itemize{
#'   \item \code{theta}: Fitted kernel lengthscale parameter(s). For TwinBKP and
#'     TwinDKP, this is an alias for the global lengthscale \code{theta_g}.
#'   \item \code{alpha_n}: Posterior Beta \eqn{\alpha} shape parameters for BKP
#'     and TwinBKP, or posterior Dirichlet concentration parameters for DKP and
#'     TwinDKP.
#'   \item \code{beta_n}: Posterior Beta \eqn{\beta} shape parameters, returned
#'     for \code{BKP} and \code{TwinBKP} objects.
#' }
#' For \code{TwinBKP} and \code{TwinDKP} objects, the returned list also includes
#' \code{theta_g}, \code{theta_l}, \code{global_kernel}, \code{local_kernel},
#' \code{global_indices}, and \code{control}.
#'
#' @keywords BKP
#'
#' @seealso \code{\link{fit_BKP}}, \code{\link{fit_DKP}},
#'   \code{\link{fit_TwinBKP}}, and \code{\link{fit_TwinDKP}} for model fitting;
#'   \code{\link{fitted}} for posterior mean extraction.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}. arXiv. <doi:10.48550/arXiv.2508.10447>.
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
#' model <- fit_BKP(X, y, m, Xbounds = Xbounds, theta = 0.3)
#'
#' # Extract posterior parameters
#' parameter(model)
#'
#' \dontrun{
#' # Larger TwinBKP example
#' n <- 200
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit TwinBKP model
#' model <- fit_TwinBKP(
#'      X, y, m,
#'      Xbounds = Xbounds,
#'      theta_g = 0.3,
#'      g = 20,
#'      twins = 1,
#'      n_threads = 1
#'    )
#'
#' # Extract posterior and kernel parameters
#' parameter(model)
#' }
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
