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
  out <- object$mean_global
  if (is.null(out)) stop("'mean_global' is missing in object. Please refit using latest fit_TwinDKP().")
  rownames(out) <- paste0("idx_", object$global_idx)
  out
}
