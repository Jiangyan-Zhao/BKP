#' @rdname fitted
#' @keywords TwinDKP
#' @export
#' @method fitted TwinDKP
fitted.TwinDKP <- function(object, ...) {
  object$alpha_n / rowSums(object$alpha_n)
}
