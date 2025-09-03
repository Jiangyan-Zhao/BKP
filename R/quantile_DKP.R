#' @rdname quantile
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
#' m <- sample(100, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit DKP model
#' model <- fit_DKP(X, Y, Xbounds = Xbounds)
#'
#' # Extract posterior quantiles
#' quantile(model)
#'
#' @export
#' @method quantile DKP

quantile.DKP <- function(x, probs = c(0.025, 0.5, 0.975), ...) {
  # arguments checking
  if (!is.numeric(probs) || any(probs < 0 | probs > 1)) {
    stop("'probs' must be a numeric vector with all values in [0, 1].")
  }

  # Extract posterior Dirichlet parameters
  alpha_n <- x$alpha_n  # n x q matrix (n obs, q classes)
  n <- nrow(alpha_n)
  q <- ncol(alpha_n)
  row_sum <- rowSums(alpha_n)

  if (length(probs) > 1) {
    # Create 3D array: obs x class x probs
    post_q_array <- array(NA, dim = c(n, q, length(probs)),
                          dimnames = list(NULL,
                                          paste0("class", 1:q),
                                          paste0(probs*100, "%")))

    # Loop over classes to compute Beta approximation quantiles
    for (j in 1:q) {
      post_q_array[, j, ] <- t(mapply(function(alpha_ij, row_sum_i) {
        qbeta(probs, alpha_ij, row_sum_i - alpha_ij)
      }, alpha_n[, j], row_sum))
    }

    return(post_q_array)
  } else {
    # Single probability: return matrix obs x class
    post_q_matrix <- matrix(NA, nrow = n, ncol = q,
                            dimnames = list(NULL, paste0("class", 1:q)))

    for (j in 1:q) {
      post_q_matrix[, j] <- mapply(function(alpha_ij, row_sum_i) {
        qbeta(probs, alpha_ij, row_sum_i - alpha_ij)
      }, alpha_n[, j], row_sum)
    }

    return(post_q_matrix)
  }
}
