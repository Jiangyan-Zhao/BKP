#' @rdname quantile
#'
#' @keywords BKP TwinBKP
#'
#' @examples
#' # -------------------------- TwinBKP ---------------------------
#' set.seed(123)
#'
#' # Define true success probability function
#' true_pi_fun <- function(x) {
#'   (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' }
#'
#' n <- 1000
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit TwinBKP model
#' model <- fit_TwinBKP(X, y, m, Xbounds = Xbounds)
#'
#' # Extract posterior quantiles
#' quantile(model)
#'
#' @export
#' @method quantile TwinBKP
quantile.TwinBKP <- function(x, probs = c(0.025, 0.5, 0.975), ...) {
  if (!is.numeric(probs) || anyNA(probs) || any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1)) {
    stop("'probs' must be a finite numeric vector with all values in [0, 1].")
  }

  alpha_n <- x$alpha_n
  beta_n <- x$beta_n
  n <- length(alpha_n)

  if (length(probs) > 1) {
    post_q <- matrix(
      stats::qbeta(
        rep(probs, each = n),
        rep(alpha_n, times = length(probs)),
        rep(beta_n, times = length(probs))
      ),
      nrow = n,
      ncol = length(probs)
    )
    colnames(post_q) <- paste0(probs * 100, "%")
  } else {
    post_q <- stats::qbeta(probs, alpha_n, beta_n)
  }

  post_q
}
