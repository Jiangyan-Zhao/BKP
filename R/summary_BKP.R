#' @name summary
#'
#' @title Summary of a Fitted BKP or DKP Model
#'
#' @description Provides a summary of a fitted Beta Kernel Process (BKP) or
#'   Dirichlet Kernel Process (DKP) model. Currently, this function acts as a
#'   wrapper for \code{\link{print.BKP}} or \code{\link{print.DKP}}, delivering
#'   a concise overview of key model characteristics and fitting results.
#'
#' @param object An object of class \code{"BKP"} (from \code{\link{fit.BKP}}) or
#'   \code{"DKP"} (from \code{\link{fit.DKP}}).
#' @param ... Additional arguments passed to the generic \code{summary} method
#'   (currently not used).
#'
#' @return Invisibly returns the input object (of class \code{"BKP"} or
#'   \code{"DKP"}). Called for side effects: prints a concise summary of the
#'   fitted model.
#'
#' @seealso \code{\link{fit.BKP}}, \code{\link{fit.DKP}},
#'   \code{\link{print.BKP}}, \code{\link{print.DKP}}.
#'
#' @keywords BKP
#'
#' @examples
#' # ============================================================== #
#' # ========================= BKP Examples ======================= #
#' # ============================================================== #
#'
#' #-------------------------- 1D Example ---------------------------
#' set.seed(123)
#'
#' # Define true success probability function
#' true_pi_fun <- function(x) {
#'   (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' }
#'
#' n <- 30
#' Xbounds <- matrix(c(-2,2), nrow=1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit BKP model
#' model1 <- fit.BKP(X, y, m, Xbounds=Xbounds)
#' summary(model1)
#'
#'
#' #-------------------------- 2D Example ---------------------------
#' set.seed(123)
#'
#' # Define 2D latent function and probability transformation
#' true_pi_fun <- function(X) {
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
#'   return(pnorm(f))  # Transform to probability
#' }
#'
#' n <- 100
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit BKP model
#' model2 <- fit.BKP(X, y, m, Xbounds=Xbounds)
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

