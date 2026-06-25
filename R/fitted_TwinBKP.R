#' @rdname fitted
#'
#' @keywords BKP TwinBKP
#'
#' @examples
#' # -------------------------- TwinBKP ---------------------------
#' set.seed(123)
#'
#' # Define true success probability function
#' true_pi_fun <- function(x) {
#'   (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' }
#'
#' n <- 1000
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit TwinBKP model
#' model <- fit_TwinBKP(X, y, m, Xbounds = Xbounds)
#'
#' # Extract fitted values
#' fitted(model)
#'
#' @export
#' @method fitted TwinBKP

fitted.TwinBKP <- function(object, ...) {
  # Posterior beta parameters
  alpha_n <- object$alpha_n
  beta_n  <- object$beta_n

  # Posterior mean
  fitted_value <- alpha_n / (alpha_n + beta_n)

  return(fitted_value)
}
