#' @rdname quantile
#' @keywords BKP
#'
#' @details
#' For \code{TwinBKP}, posterior quantiles are computed by first calling
#' \code{predict.TwinBKP()} at \code{Xnew}, then applying \code{qbeta()} to
#' returned \code{alpha_n} and \code{beta_n}. If \code{Xnew = NULL}, the method
#' generates \code{n_grid} points in \code{Xbounds} via Latin hypercube sampling
#' and then performs prediction.
#'
#' @param Xnew Optional prediction points for quantile evaluation.
#' @param n_grid Number of generated points when \code{Xnew = NULL}. Default is
#'   \code{1000}.
#' @param l_nums Optional local neighbor count passed to
#'   \code{predict.TwinBKP()}.
#' @param v_nums Optional validation size passed to
#'   \code{predict.TwinBKP()}.
#'
#' @export
#' @method quantile TwinBKP
quantile.TwinBKP <- function(x, probs = c(0.025, 0.5, 0.975),
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

  pred <- predict.TwinBKP(
    object = x,
    Xnew = Xnew,
    l_nums = l_nums,
    v_nums = v_nums,
    ...
  )

  alpha_n <- pred$alpha_n
  beta_n  <- pred$beta_n

  if (length(probs) > 1) {
    post_q <- t(mapply(function(a, b) qbeta(probs, a, b), alpha_n, beta_n))
    colnames(post_q) <- paste0(probs * 100, "%")
  } else {
    post_q <- mapply(function(a, b) qbeta(probs, a, b), alpha_n, beta_n)
  }

  post_q
}
