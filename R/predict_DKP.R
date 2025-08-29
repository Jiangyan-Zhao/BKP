#' @rdname predict
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
#'
#' # Prediction on training data
#' predict(model1)
#'
#' # Prediction on new data
#'
#' # Prediction
#' Xnew = matrix(seq(-2, 2, length = 10), ncol=1) #new data points
#' predict(model1, Xnew)
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
#'
#' # Prediction on training data
#' predict(model2)
#'
#' # Prediction on new data
#' x1 <- seq(Xbounds[1,1], Xbounds[1,2], length.out = 10)
#' x2 <- seq(Xbounds[2,1], Xbounds[2,2], length.out = 10)
#' Xnew <- expand.grid(x1 = x1, x2 = x2)
#' predict(model2, Xnew)
#'
#' @export
#' @method predict DKP

predict.DKP <- function(object, Xnew = NULL, CI_level = 0.95, ...)
{
  # Extract components
  X       <- object$X
  Xnorm   <- object$Xnorm
  Y       <- object$Y
  theta   <- object$theta_opt
  kernel  <- object$kernel
  prior   <- object$prior
  r0      <- object$r0
  p0      <- object$p0
  Xbounds <- object$Xbounds
  d       <- ncol(Xnorm)

  if(!is.null(Xnew)){
    # Ensure Xnew is a matrix and matches input dimension
    if (is.null(nrow(Xnew))) {
      Xnew <- matrix(Xnew, nrow = 1)
    }
    Xnew <- as.matrix(Xnew)
    if (ncol(Xnew) != d) {
      stop("The number of columns in 'Xnew' must match the original input dimension.")
    }

    # Normalize Xnew to [0,1]^d
    Xnew_norm <- sweep(Xnew, 2, Xbounds[, 1], "-")
    Xnew_norm <- sweep(Xnew_norm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

    # Compute kernel matrix
    K <- kernel_matrix(Xnew_norm, Xnorm, theta = theta, kernel = kernel) # [m × n]
  }else{
    # Compute kernel matrix
    K <- kernel_matrix(Xnorm, theta = theta, kernel = kernel) # [n × n]
  }

  # get the prior parameters: alpha0(x) and beta0(x)
  alpha0 <- get_prior_dkp(prior = prior, r0 = r0, p0 = p0, Y = Y, K = K)

  # Posterior parameters
  alpha_n <- pmax(alpha0 + as.matrix(K %*% Y), 1e-6) # [n × q]

  # Predictive quantities
  row_sum <- rowSums(alpha_n)
  pi_mean <- alpha_n / row_sum
  pi_var  <- alpha_n * (row_sum - alpha_n) / (row_sum^2 * (row_sum + 1))
  pi_lower <- qbeta((1 - CI_level) / 2, alpha_n, row_sum - alpha_n)
  pi_upper <- qbeta((1 + CI_level) / 2, alpha_n, row_sum - alpha_n)

  # Return structured output
  prediction <- list(
    X        = X,
    Xnew     = Xnew,
    mean     = pi_mean,
    variance = pi_var,
    lower    = pi_lower,
    upper    = pi_upper,
    CI_level  = CI_level
  )

  # Posterior classification label (only for classification data)
  if (all(rowSums(Y) == 1)) {
    prediction$class <- max.col(pi_mean)
  }

  class(prediction) <- "predict.DKP"
  return(prediction)
}
