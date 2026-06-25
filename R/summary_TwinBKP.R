#' @rdname summary
#'
#' @keywords BKP TwinBKP
#'
#' @examples
#' # ============================================================== #
#' # ======================= TwinBKP Examples ===================== #
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
#' n <- 1000
#' Xbounds <- matrix(c(-2,2), nrow=1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit TwinBKP model
#' model1 <- fit_TwinBKP(X, y, m, Xbounds=Xbounds)
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
#' n <- 1000
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit TwinBKP model
#' model2 <- fit_TwinBKP(X, y, m, Xbounds=Xbounds)
#' summary(model2)
#'
#' @export
#' @method summary TwinBKP
summary.TwinBKP <- function(object, ...) {

  n_obs <- nrow(object$X)
  d <- ncol(object$X)

  post_mean <- as.vector(object$alpha_n / (object$alpha_n + object$beta_n))
  post_var <- as.vector(
    post_mean * (1 - post_mean) / (object$alpha_n + object$beta_n + 1)
  )

  control <- object$control

  res <- list(
    n_obs = n_obs,
    input_dim = d,

    kernel = object$kernel,
    global_kernel = object$global_kernel,
    local_kernel = object$local_kernel,

    isotropic = object$isotropic,

    theta_opt = object$theta_opt,
    theta_g = object$theta_g,
    theta_l = object$theta_l,

    loss = object$loss,
    loss_min = object$loss_min,

    prior = object$prior,
    r0 = object$r0,
    p0 = object$p0,

    global_size = length(object$global_indices),
    global_target = control$g_target,
    local_size = control$l,
    twins = control$twins,

    post_mean = post_mean,
    post_var = post_var
  )

  class(res) <- "summary_TwinBKP"
  res
}
