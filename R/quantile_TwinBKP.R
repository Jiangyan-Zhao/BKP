#' @rdname quantile
#' @keywords BKP TwinBKP
#' @export
#' @method quantile TwinBKP
quantile.TwinBKP <- function(x, probs = c(0.025, 0.5, 0.975), ...) {
  # TwinBKP has the same Beta posterior-parameter representation as BKP.
  quantile.BKP(x, probs = probs, ...)
}
