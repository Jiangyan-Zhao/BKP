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
#' @examples
#' \dontrun{
#' ### 1D
#' set.seed(123)
#' n <- 100
#' Xbounds <- matrix(c(-2,2), nrow=1)
#' x <- seq(-2, 2, length = n)
#' true_pi <- (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#' df <- data.frame(x = x, y = y, m = m)
#' xx = matrix(seq(-2, 2, length = 100), ncol=1) #new data points
#' model <- fit.BKP(df, Xbounds=Xbounds)
#' print(model)
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
#' library(tgp)
#' x <- lhs(n = n, rect = Xbounds)
#' true_pi <- pnorm(f(x))
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#' df <- data.frame(x = x, y = y, m = m)
#' xx1 <- seq(Xbounds[1,1], Xbounds[1,2], length.out = 100)
#' xx2 <- seq(Xbounds[2,1], Xbounds[2,2], length.out = 100)
#' xx <- expand.grid(xx1 = xx1, xx2 = xx2)
#' model <- fit.BKP(df)
#' print(model)
#' }
#'
#' @export
#' @method print BKP

print.BKP <- function(x, ...) {
  if (!inherits(x, "BKP")) {
    stop("The input is not of class 'BKP'. Please provide a model fitted with 'fit.BKP()'.")
  }

  n <- nrow(x$X)
  d <- ncol(x$X)
  theta <- x$bestTheta
  kernel <- x$kernel
  loss <- x$loss
  minLoss <- x$minLoss
  prior <- x$prior
  r0 <- x$r0
  p0 <- x$p0

  cat("--------------------------------------------------\n")
  cat("           Beta Kernel Process (BKP) Model        \n")
  cat("--------------------------------------------------\n")
  cat(sprintf("Number of observations (n):  %d\n", n))
  cat(sprintf("Input dimensionality (d):    %d\n", d))
  cat(sprintf("Kernel type:                 %s\n", kernel))
  cat(sprintf("Loss function used:          %s\n", loss))
  cat(sprintf("Optimized kernel parameters: %s\n",
              paste(sprintf("%.4f", theta), collapse = ", ")))
  cat(sprintf("Minimum achieved loss:       %.5f\n", minLoss))
  cat("\n")

  cat("Prior specification:\n")
  cat(sprintf("  Type:    %s\n", prior))
  if (prior == "adaptive") {
    cat("  Note:    Data-adaptive informative prior used.\n")
    cat(sprintf("  r0:      %.3f\n", r0))
  } else if (prior == "fixed") {
    cat("  Note:    Fixed informative prior shared across locations.\n")
    cat(sprintf("  r0:      %.3f\n", r0))
    cat(sprintf("  p0:      %.3f\n", p0))
  } else if (prior == "noninformative") {
    cat("  Note:    Noninformative prior: Beta(1,1).\n")
  }

  cat("--------------------------------------------------\n")
}

