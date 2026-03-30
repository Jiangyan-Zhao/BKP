 #' @rdname fitted
#' @keywords BKP
#'
#' @details
#' For \code{TwinBKP}, fitted values are computed on the global support set
#' selected during \code{fit_TwinBKP()}. Specifically, this method computes
#' posterior Beta parameters on \code{object$Xnorm_global} and returns the
#' posterior mean \eqn{\alpha_n / (\alpha_n + \beta_n)} for those support
#' points.
#'
#' This is intentionally different from \code{BKP}, where fitted values are
#' naturally available on all training points. In TwinBKP, full-sample behavior
#' is prediction-oriented and should be obtained via \code{predict.TwinBKP()}.
#'
#' @export
#' @method fitted TwinBKP
fitted.TwinBKP <- function(object, ...) {
  gp <- .twinbkp_global_posterior(object)
  out <- gp$mean
  names(out) <- paste0("idx_", object$global_idx)
  out
}
