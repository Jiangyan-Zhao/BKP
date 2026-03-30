#' @rdname fitted
#' @keywords DKP
#'
#' @details
#' For \code{TwinDKP}, fitted values are reported on the global support points
#' selected during \code{fit_TwinDKP()}. The method computes posterior mean
#' class probabilities on \code{object$Xnorm_global}.
#'
#' @export
#' @method fitted TwinDKP
fitted.TwinDKP <- function(object, ...) {
  gp <- .twindkp_global_posterior(object)
  out <- gp$mean
  rownames(out) <- paste0("idx_", object$global_idx)
  out
}
