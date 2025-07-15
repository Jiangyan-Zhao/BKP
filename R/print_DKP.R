#' @name print
#'
#' @title Print Summary of a DKP Model
#'
#' @description
#' Provides a concise summary of a fitted Beta Kernel Process (DKP) model object.
#' The printed output includes key model characteristics such as the number of observations,
#' the dimensionality of the input, the kernel and loss type used, the optimized kernel
#' parameters (`bestTheta`), and the minimum loss (e.g., Brier score) achieved.
#'
#' @param x An object of class \code{"DKP"}, typically returned by \code{\link{fit.DKP}}.
#' @param ... Additional arguments passed to the generic \code{print} method (not used here).
#'
#' @details
#' When a \code{DKP} object is printed (e.g., by calling \code{print(model)} or typing the object name
#' in the console), this method presents a readable summary of the fitted model's essential properties.
#'
#' @seealso
#' \code{\link{fit.DKP}} for model fitting.
#' \code{\link{summary.DKP}} for a more detailed report, if implemented.
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
#' print(DKPmodel)
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
#' print(DKPmodel)
#' }
#'
#' @export
#' @method print DKP

print.DKP <- function(x, ...) {
  if (!inherits(x, "DKP")) {
    stop("The input is not of class 'DKP'. Please provide a model fitted with 'fit.DKP()'.")
  }

  n <- nrow(x$X)
  d <- ncol(x$X)
  q <- ncol(x$Y)
  theta <- x$bestTheta
  kernel <- x$kernel
  loss <- x$loss
  minLoss <- x$minLoss
  prior <- x$prior
  r0 <- x$r0
  p0 <- x$p0

  cat("--------------------------------------------------\n")
  cat("           Dirichlet Kernel Process (DKP) Model        \n")
  cat("--------------------------------------------------\n")
  cat(sprintf("Number of observations (n):  %d\n", n))
  cat(sprintf("Input dimensionality (d):    %d\n", d))
  cat(sprintf("Output dimensionality (q):    %d\n", q))
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
    cat("  Note:    Noninformative prior: Dirichlet(1,...,1).\n")
  }

  cat("--------------------------------------------------\n")
}

