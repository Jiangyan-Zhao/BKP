#' @name summary
#'
#' @title Summary of a DKP Model
#'
#' @description
#' Provides a summary of a fitted Beta Kernel Process (DKP) model.
#' Currently, this function serves as a wrapper for \code{\link{print.DKP}},
#' offering a concise overview of the model's key components and fitting results.
#'
#' @param object An object of class \code{"DKP"}, typically returned by \code{\link{fit.DKP}}.
#' @param ... Additional arguments passed to the generic \code{summary} method (not used here).
#'
#' @details
#' The \code{summary.DKP} method displays essential model details such as the number of
#' observations, input dimensionality, kernel type, optimized kernel parameters,
#' and the corresponding loss value (e.g., Brier score). Future versions may expand this
#' method to include more diagnostics, confidence metrics, or model evaluation summaries.
#'
#' @seealso
#' \code{\link{print.DKP}} for basic print output.
#' \code{\link{fit.DKP}} to fit a DKP model.
#'
#' @author Jiangyan Zhao, Kunhai Qing, Jin Xu
#'
#' @keywords DKP
#'
#' @examples
#' \dontrun{
#' ### 1D
#' set.seed(123)
#' n <- 30
#' Xbounds <- matrix(c(-2,2), nrow=1)
#' x <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' true_pi <- matrix(c(true_pi/2,true_pi/2,1-true_pi),nrow = n, byrow = F)
#' m <- sample(100, n, replace = TRUE)
#' Y <- matrix(0, nrow = n, ncol = 3)
#' for (i in 1:n) {
#'   Y[i, ] <- rmultinom(n=1, size=m[i], prob=true_pi[i, ])
#' }
#' DKPmodel <- fit.DKP(p=3,X=x,Y=Y,Xbounds = Xbounds,prior = "noninformative",kernel = "gaussian",loss = "brier")
#' summary(DKPmodel)
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
#' true_pi <- pnorm(f(x))
#' true_pi <- matrix(c(true_pi/2,true_pi/2,1-true_pi),nrow = n, byrow = F)
#' m <- sample(100, n, replace = TRUE)
#' Y <- matrix(0, nrow = n, ncol = 3)
#' for (i in 1:n) {
#'   Y[i, ] <- rmultinom(n=1, size=m[i], prob=true_pi[i, ])
#' }
#' DKPmodel <- fit.DKP(p=3,X=x,Y=Y,Xbounds = Xbounds,prior = "noninformative",kernel = "gaussian",loss = "brier")
#' summary(DKPmodel)
#' }
#' @export
#' @method summary DKP

summary.DKP <- function(object, ...) {
  if (!inherits(object, "DKP")) {
    stop("The input is not of class 'DKP'. Please provide a fitted DKP model.")
  }

  # Delegate to print.DKP for now
  print(object, ...)
}

