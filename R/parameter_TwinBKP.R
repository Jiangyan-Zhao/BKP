#' @rdname parameter
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
#' n <- 30
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit TwinBKP model
#' model <- fit_TwinBKP(X, y, m, Xbounds = Xbounds)
#'
#' # Extract posterior and kernel parameters
#' parameter(model)
#'
#' @export
#' @method parameter TwinBKP
parameter.TwinBKP <- function(object, ...) {
  list(
    theta = object$theta_opt,
    theta_g = object$theta_g,
    theta_l = object$theta_l,
    global_kernel = object$global_kernel,
    local_kernel = object$local_kernel,
    alpha_n = object$alpha_n,
    beta_n = object$beta_n,
    global_indices = object$global_indices,
    control = object$control
  )
}
