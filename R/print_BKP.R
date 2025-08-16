#' @name print
#'
#' @title Print Summary of a Fitted BKP or DKP Model
#'
#' @description Displays a concise summary of a fitted BKP or DKP model. The
#'   output includes key characteristics such as sample size, input
#'   dimensionality, kernel type, loss function, optimized kernel
#'   hyperparameters, and minimum loss.
#'
#' @param x An object of class \code{"BKP"} (from \code{\link{fit.BKP}}) or
#'   \code{"DKP"} (from \code{\link{fit.DKP}}).
#' @param ... Additional arguments passed to the generic \code{print} method
#'   (currently not used).
#'
#' @return Invisibly returns the input object (of class \code{"BKP"} or
#'   \code{"DKP"}). The function is called for its side effect of printing a
#'   summary to the console.
#'
#' @seealso \code{\link{fit.BKP}}, \code{\link{fit.DKP}},
#'   \code{\link{summary.BKP}}, \code{\link{summary.DKP}}.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv.
#'   https://doi.org/10.48550/arXiv.2508.10447.
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
#' print(model1)
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
#' print(model2)
#'
#' @export
#' @method print BKP

print.BKP <- function(x, ...) {
  if (!inherits(x, "BKP")) {
    stop("The input is not of class 'BKP'. Please provide a model fitted with 'fit.BKP()'.")
  }

  n <- nrow(x$X)
  d <- ncol(x$X)
  theta <- x$theta_opt
  kernel <- x$kernel
  loss <- x$loss
  loss_min <- x$loss_min
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
  if (!is.na(loss_min)) {
    cat(sprintf("Minimum achieved loss:       %.5f\n", loss_min))
    cat("Kernel parameters were obtained by optimization.\n")
  } else {
    cat("Note: Kernel parameters were user-specified (no optimization).\n")
  }
  cat("\n")

  cat("Prior specification:\n")
  if (prior == "adaptive") {
    cat("  Data-adaptive informative prior used.\n")
    cat(sprintf("  r0:      %.3f\n", r0))
  } else if (prior == "fixed") {
    cat("  Fixed informative prior shared across locations.\n")
    cat(sprintf("  r0:      %.3f\n", r0))
    cat(sprintf("  p0:      %.3f\n", p0))
  } else if (prior == "noninformative") {
    cat("  Noninformative prior: Beta(1,1).\n")
  }

  cat("--------------------------------------------------\n")
}

