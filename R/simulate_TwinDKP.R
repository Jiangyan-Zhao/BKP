#' @rdname simulate
#' @keywords DKP
#'
#' @details
#' For \code{TwinDKP}, posterior simulation is prediction-driven:
#' \enumerate{
#'   \item Build/receive \code{Xnew}.
#'   \item Call \code{predict.TwinDKP()} to obtain \code{alpha_n}.
#'   \item Draw \code{nsim} samples from class-wise Dirichlet posteriors.
#' }
#'
#' If \code{Xnew = NULL}, \code{n_grid} points are generated in \code{Xbounds}
#' (default \code{1000}) and used for simulation.
#'
#' @param n_grid Number of generated points when \code{Xnew = NULL}. Default is
#'   \code{1000}.
#' @param l_nums Optional local neighbor count passed to
#'   \code{predict.TwinDKP()}.
#' @param v_nums Optional validation size passed to
#'   \code{predict.TwinDKP()}.
#'
#' @export
#' @method simulate TwinDKP
simulate.TwinDKP <- function(object, nsim = 1, seed = NULL,
                             Xnew = NULL, n_grid = 1000,
                             l_nums = NULL, v_nums = NULL, ...) {
  if (!is.numeric(nsim) || length(nsim) != 1 || nsim <= 0 || nsim != as.integer(nsim)) {
    stop("`nsim` must be a positive integer.")
  }
  nsim <- as.integer(nsim)

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1 || seed != as.integer(seed))) {
    stop("`seed` must be a single integer or NULL.")
  }
  if (!is.null(seed)) set.seed(seed)

  if (is.null(Xnew)) {
    if (!is.numeric(n_grid) || length(n_grid) != 1 || n_grid <= 0) {
      stop("'n_grid' must be a positive integer when 'Xnew' is NULL.")
    }
    Xnew <- tgp::lhs(as.integer(n_grid), rect = object$Xbounds)
  }

  pred <- predict.TwinDKP(
    object = object,
    Xnew = Xnew,
    l_nums = l_nums,
    v_nums = v_nums,
    ...
  )

  alpha_n <- as.matrix(pred$alpha_n)
  n_new <- nrow(alpha_n)
  q <- ncol(alpha_n)

  samples <- array(0, dim = c(n_new, q, nsim))
  for (i in seq_len(n_new)) {
    samples[i, , ] <- t(dirmult::rdirichlet(n = nsim, alpha = alpha_n[i, ]))
  }
  dimnames(samples) <- list(
    paste0("x", seq_len(n_new)),
    colnames(alpha_n),
    paste0("sim", seq_len(nsim))
  )

  class_pred <- NULL
  if (!is.null(pred$class)) {
    class_pred <- matrix(NA_integer_, nrow = n_new, ncol = nsim)
    for (i in seq_len(nsim)) class_pred[, i] <- max.col(samples[, , i])
    colnames(class_pred) <- paste0("sim", seq_len(nsim))
  }

  res <- list(
    samples = samples,
    mean = as.matrix(pred$mean),
    class = class_pred,
    X = object$X,
    Xnew = Xnew
  )
  class(res) <- "simulate_TwinDKP"
  res
}
