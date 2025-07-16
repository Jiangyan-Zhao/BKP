#' @name summary
#'
#' @title Summary of a BKP Model
#'
#' @description Provides a summary of a fitted Beta Kernel Process (BKP) model.
#'   Currently, this function serves as a wrapper for \code{\link{print.BKP}},
#'   offering a concise overview of the model's key components and fitting
#'   results.
#'
#' @param object An object of class \code{"BKP"}, typically returned by
#'   \code{\link{fit.BKP}}.
#' @param ... Additional arguments passed to the generic \code{summary} method
#'   (not used here).
#'
#' @details The \code{summary.BKP} method displays essential model details such
#'   as the number of observations, input dimensionality, kernel type, optimized
#'   kernel parameters, and the corresponding loss value (e.g., Brier score).
#'   Future versions may expand this method to include more diagnostics,
#'   confidence metrics, or model evaluation summaries.
#'
#' @seealso \code{\link{print.BKP}} for basic print output.
#'   \code{\link{fit.BKP}} to fit a BKP model.
#'
#' @author Jiangyan Zhao, Kunhai Qing, Jin Xu
#'
#' @keywords BKP
#'
#' @examples
#' ### 1D
#' set.seed(123)
#' n <- 30
#' Xbounds <- matrix(c(-2,2), nrow=1)
#' x <- tgp::lhs(n = n, rect = Xbounds)
#' m <- sample(100, n, replace = TRUE)
#' true_pi <- (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' y <- rbinom(n, size = m, prob = true_pi)
#' model1 <- fit.BKP(x, y, m, Xbounds=Xbounds)
#' summary(model1)
#'
#' ### 2D
#' set.seed(123)
#' n <- 100
#' f <- function(X) {
#'   if(is.null(nrow(X))) X <- matrix(X, nrow=1)
#'   m <- 8.6928
#'   s <- 2.4269
#'   x1 <- 4*X[,1]- 2
#'   x2 <- 4*X[,2]- 2
#'   a <- 1 + (x1 + x2 + 1)^2 *
#'     (19- 14*x1 + 3*x1^2- 14*x2 + 6*x1*x2 + 3*x2^2)
#'   b <- 30 + (2*x1- 3*x2)^2 *
#'     (18- 32*x1 + 12*x1^2 + 48*x2- 36*x1*x2 + 27*x2^2)
#'   f <- log(a*b)
#'   f <- (f- m)/s
#'   return(f) }
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' x <- tgp::lhs(n = n, rect = Xbounds)
#' m <- sample(100, n, replace = TRUE)
#' true_pi <- pnorm(f(x))
#' y <- rbinom(n, size = m, prob = true_pi)
#' model2 <- fit.BKP(x, y, m)
#' summary(model2)
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

