#' @name summary
#'
#' @title Summary of a BKP Model
#'
#' @description
#' Provides a summary of a fitted Beta Kernel Process (BKP) model.
#' Currently, this function serves as a wrapper for \code{\link{print.BKP}},
#' offering a concise overview of the model's key components and fitting results.
#'
#' @param object An object of class \code{"BKP"}, typically returned by \code{\link{fit.BKP}}.
#' @param ... Additional arguments passed to the generic \code{summary} method (not used here).
#'
#' @details
#' The \code{summary.BKP} method displays essential model details such as the number of
#' observations, input dimensionality, kernel type, optimized kernel parameters,
#' and the corresponding loss value (e.g., Brier score). Future versions may expand this
#' method to include more diagnostics, confidence metrics, or model evaluation summaries.
#'
#' @seealso
#' \code{\link{print.BKP}} for basic print output.
#' \code{\link{fit.BKP}} to fit a BKP model.
#'
#' @author Jiangyan Zhao, Kunhai Qing, Jin Xu
#'
#' @keywords BKP
#'
#' @export
#' @method summary BKP

summary.BKP <- function(object, ...) {
  if (!inherits(object, "BKP")) {
    stop("The input is not of class 'BKP'. Please provide a fitted BKP model.")
  }

  # Delegate to print.BKP for now
  print(object, ...)
}

