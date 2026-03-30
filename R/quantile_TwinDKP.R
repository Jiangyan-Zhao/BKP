#' @rdname quantile
#' @keywords DKP
#'
#' @details
#' For \code{TwinDKP}, posterior quantiles are prediction-driven.
#' The method first calls \code{predict.TwinDKP()} at \code{Xnew} and then
#' computes class-wise marginal quantiles via Beta approximation:
#' for class \eqn{j}, use \code{qbeta(probs, alpha_{ij}, sum(alpha_i)-alpha_{ij})}.
#'
#' If \code{Xnew = NULL}, \code{n_grid} points are generated in
#' \code{x$Xbounds} by Latin hypercube sampling.
#'
#' @param Xnew Optional prediction points for quantile evaluation.
#' @param n_grid Number of generated points when \code{Xnew = NULL}. Default is
#'   \code{1000}.
#' @param l_nums Optional local neighbor count passed to
#'   \code{predict.TwinDKP()}.
#' @param v_nums Optional validation size passed to
#'   \code{predict.TwinDKP()}.
#'
#' @export
#' @method quantile TwinDKP
quantile.TwinDKP <- function(x, probs = c(0.025, 0.5, 0.975),
                             Xnew = NULL, n_grid = 1000,
                             l_nums = NULL, v_nums = NULL, ...) {
  if (!is.numeric(probs) || any(probs < 0 | probs > 1)) {
    stop("'probs' must be a numeric vector with all values in [0, 1].")
  }

  if (is.null(Xnew)) {
    if (!is.numeric(n_grid) || length(n_grid) != 1 || n_grid <= 0) {
      stop("'n_grid' must be a positive integer when 'Xnew' is NULL.")
    }
    Xnew <- tgp::lhs(as.integer(n_grid), rect = x$Xbounds)
  }

  pred <- predict.TwinDKP(
    object = x,
    Xnew = Xnew,
    l_nums = l_nums,
    v_nums = v_nums,
    ...
  )

  alpha_n <- as.matrix(pred$alpha_n)
  n <- nrow(alpha_n)
  q <- ncol(alpha_n)
  row_sum <- rowSums(alpha_n)

  if (length(probs) > 1) {
    post_q_array <- array(NA, dim = c(n, q, length(probs)),
                          dimnames = list(NULL, colnames(alpha_n), paste0(probs * 100, "%")))

    for (j in seq_len(q)) {
      post_q_array[, j, ] <- t(mapply(function(alpha_ij, row_sum_i) {
        qbeta(probs, alpha_ij, pmax(row_sum_i - alpha_ij, 1e-10))
      }, alpha_n[, j], row_sum))
    }
    return(post_q_array)
  }

  post_q_matrix <- matrix(NA, nrow = n, ncol = q,
                          dimnames = list(NULL, colnames(alpha_n)))
  for (j in seq_len(q)) {
    post_q_matrix[, j] <- mapply(function(alpha_ij, row_sum_i) {
      qbeta(probs, alpha_ij, pmax(row_sum_i - alpha_ij, 1e-10))
    }, alpha_n[, j], row_sum)
  }
  post_q_matrix
}
