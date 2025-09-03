#' @rdname summary
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
#'   p1 <- 1/(1+exp(-3*X))
#'   p2 <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
#'   return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
#' }
#'
#' n <- 30
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(150, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit DKP model
#' model1 <- fit_DKP(X, Y, Xbounds = Xbounds)
#' summary(model1)
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
#'   f <- (log(a*b)- m)/s
#'   p1 <- pnorm(f) # Transform to probability
#'   p2 <- sin(pi * X[,1]) * sin(pi * X[,2])
#'   return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
#' }
#'
#' n <- 100
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(150, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit DKP model
#' model2 <- fit_DKP(X, Y, Xbounds = Xbounds)
#' summary(model2)
#'
#' @export
#' @method summary DKP

summary.DKP <- function(object, ...) {
  # Extract info
  n_obs <- nrow(object$X)
  d     <- ncol(object$X)
  q     <- ncol(object$Y)

  # Posterior Dirichlet parameters
  alpha_n <- object$alpha_n
  row_sum <- rowSums(alpha_n)

  # Posterior mean and variance for each category at each training point
  post_mean <- alpha_n / row_sum
  post_var <- (alpha_n * (row_sum - alpha_n)) / (row_sum^2 * (row_sum + 1))

  res <- list(
    n_obs     = n_obs,
    input_dim = d,
    n_class   = q,
    kernel    = object$kernel,
    theta_opt = object$theta_opt,
    loss      = object$loss,
    loss_min  = object$loss_min,
    prior     = object$prior,
    r0        = object$r0,
    p0        = object$p0,
    post_mean = post_mean,
    post_var  = post_var
  )
  class(res) <- "summary_DKP"
  return(res)
}
