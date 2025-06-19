#' @name print
#'
#' @title Print Summary of a BKP Model
#'
#' @description
#' Provides a concise summary of a fitted Beta Kernel Process (BKP) model object.
#' The printed output includes key model characteristics such as the number of observations,
#' the dimensionality of the input, the kernel and loss type used, the optimized kernel
#' parameters (`bestTheta`), and the minimum loss (e.g., Brier score) achieved.
#'
#' @param x An object of class \code{"BKP"}, typically returned by \code{\link{fit.BKP}}.
#' @param ... Additional arguments passed to the generic \code{print} method (not used here).
#'
#' @details
#' When a \code{BKP} object is printed (e.g., by calling \code{print(model)} or typing the object name
#' in the console), this method presents a readable summary of the fitted model's essential properties.
#'
#' @seealso
#' \code{\link{fit.BKP}} for model fitting.
#' \code{\link{summary.BKP}} for a more detailed report, if implemented.
#'
#' @author Jiangyan Zhao, Kunhai Qing, Jin Xu
#'
#' @keywords BKP
#'
#' @export
#' @method print BKP

print.BKP <- function(x, ...) {
  if (!inherits(x, "BKP")) {
    stop("The input is not of class 'BKP'. Please provide a model fitted with 'fit.BKP()'.")
  }

  BKPmodel <- x

  cat("--------------------------------------------------\n")
  cat("         Beta Kernel Process (BKP) Summary        \n")
  cat("--------------------------------------------------\n\n")

  cat("Number of Observations (n): ", nrow(BKPmodel$X), "\n")
  cat("Input Dimensionality (d):   ", ncol(BKPmodel$X), "\n")
  cat("Kernel Type:                ", BKPmodel$kernel, "\n")
  cat("Loss Type:                  ", BKPmodel$loss, "\n")
  cat("Optimized Parameters (θ):   ", paste(round(BKPmodel$bestTheta, 4), collapse = ", "), "\n")
  cat("Minimum Loss Value:         ", round(BKPmodel$minLoss, 4), "\n\n")
}
