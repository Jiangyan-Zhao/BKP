#' @rdname fitted
#' @keywords DKP
#'
#' @examples
#' # -------------------------- DKP ---------------------------
#' #' set.seed(123)
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
#' model <- fit_DKP(X, Y, Xbounds = Xbounds)
#'
#' # Extract fitted values
#' fitted(model)
#'
#' @export
#' @method fitted DKP

fitted.DKP <- function(object, ...) {
  # assume n x q matrix (q classes)
  alpha_n <- object$alpha_n
  q <- ncol(alpha_n)

  # Posterior mean
  fitted_value <- alpha_n / rowSums(alpha_n)
  colnames(fitted_value) <- paste0("class", seq_len(q))

  return(fitted_value)
}
