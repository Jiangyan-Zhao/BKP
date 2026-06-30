#' @rdname quantile
#' @keywords TwinDKP
#' @export
#' @method quantile TwinDKP
quantile.TwinDKP <- function(x, probs = c(0.025, 0.5, 0.975), ...) {
  # TwinDKP has the same Dirichlet posterior-parameter representation as DKP.
  quantile.DKP(x, probs = probs, ...)
}
