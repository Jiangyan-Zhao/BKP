#' @rdname fitted
#' @keywords TwinBKP
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
