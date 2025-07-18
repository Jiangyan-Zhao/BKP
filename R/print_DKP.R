#' @rdname print
#'
#' @keywords DKP
#'
#' @examples
#' # ============================================================== #
#' # ========================= DKP Examples ======================= #
#' # ============================================================== #
#'
#' #-------------------------- 1D Example ---------------------------
#' set.seed(123)
#'
#' # Define true class probability function (3-class)
#' true_pi_fun <- function(X) {
#'   p <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
#'   return(matrix(c(p/2, p/2, 1 - p), nrow = length(p)))
#' }
#'
#' n <- 30
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit DKP model
#' model1 <- fit.DKP(X, Y, Xbounds = Xbounds)
#' print(model1)
#'
#'
#' #-------------------------- 2D Example ---------------------------
#' set.seed(123)
#'
#' # Define latent function and transform to 3-class probabilities
#' true_pi_fun <- function(X) {
#'   if (is.null(nrow(X))) X <- matrix(X, nrow = 1)
#'   m <- 8.6928; s <- 2.4269
#'   x1 <- 4 * X[,1] - 2
#'   x2 <- 4 * X[,2] - 2
#'   a <- 1 + (x1 + x2 + 1)^2 *
#'     (19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2)
#'   b <- 30 + (2*x1 - 3*x2)^2 *
#'     (18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2)
#'   f <- (log(a * b) - m) / s
#'   p <- pnorm(f)
#'   return(matrix(c(p/2, p/2, 1 - p), nrow = length(p)))
#' }
#'
#' n <- 100
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit DKP model
#' model2 <- fit.DKP(X, Y, Xbounds = Xbounds)
#' print(model2)
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
  theta <- x$theta_opt
  kernel <- x$kernel
  loss <- x$loss
  loss_min <- x$loss_min
  prior <- x$prior
  r0 <- x$r0
  p0 <- x$p0

  cat("--------------------------------------------------\n")
  cat("           Dirichlet Kernel Process (DKP) Model        \n")
  cat("--------------------------------------------------\n")
  cat(sprintf("Number of observations (n):  %d\n", n))
  cat(sprintf("Input dimensionality (d):    %d\n", d))
  cat(sprintf("Output dimensionality (q):   %d\n", q))
  cat(sprintf("Kernel type:                 %s\n", kernel))
  cat(sprintf("Loss function used:          %s\n", loss))
  cat(sprintf("Optimized kernel parameters: %s\n",
              paste(sprintf("%.4f", theta), collapse = ", ")))
  cat(sprintf("Minimum achieved loss:       %.5f\n", loss_min))
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
    cat("  Noninformative prior: Dirichlet(1,...,1).\n")
  }

  cat("--------------------------------------------------\n")
}

