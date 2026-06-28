#' @rdname quantile
#' @keywords TwinDKP
#' @export
#' @method quantile TwinDKP
quantile.TwinDKP <- function(x, probs = c(0.025, 0.5, 0.975), ...) {
  quantile.DKP(x, probs = probs, ...)
}
