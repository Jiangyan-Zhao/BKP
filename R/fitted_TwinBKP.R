#' @rdname fitted
#'
#' @keywords BKP TwinBKP
#'
#' @export
#' @method fitted TwinBKP
fitted.TwinBKP <- function(object, ...) {
  object$alpha_n / (object$alpha_n + object$beta_n)
}
