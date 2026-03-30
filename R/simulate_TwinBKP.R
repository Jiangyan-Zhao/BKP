#' @rdname simulate
#' @keywords BKP
#'
#' @details
#' For \code{TwinBKP}, posterior simulation is prediction-driven:
#' \enumerate{
#'   \item Build/receive \code{Xnew}.
#'   \item Call \code{predict.TwinBKP()} to obtain \code{alpha_n, beta_n}.
#'   \item Draw \code{nsim} samples from \code{Beta(alpha_n, beta_n)} at each
#'   point.
#' }
#'
#' If \code{Xnew = NULL}, \code{n_grid} points are generated in \code{Xbounds}
#' (default \code{1000}) and used for simulation.
#'
#' @param n_grid Number of generated points when \code{Xnew = NULL}. Default is
#'   \code{1000}.
#' @param l_nums Optional local neighbor count passed to
#'   \code{predict.TwinBKP()}.
#' @param v_nums Optional validation size passed to
#'   \code{predict.TwinBKP()}.
#'
#' @export
#' @method simulate TwinBKP
simulate.TwinBKP <- function(object, nsim = 1, seed = NULL,
                             Xnew = NULL, n_grid = 1000,
                             l_nums = NULL, v_nums = NULL,
                             threshold = NULL, ...) {
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

  pred <- predict.TwinBKP(
    object = object,
    Xnew = Xnew,
    l_nums = l_nums,
    v_nums = v_nums,
    ...
  )

  alpha_n <- pred$alpha_n
  beta_n  <- pred$beta_n
  n_new <- nrow(Xnew)

  samples <- matrix(
    rbeta(n_new * nsim,
          shape1 = rep(alpha_n, nsim),
          shape2 = rep(beta_n, nsim)),
    nrow = n_new,
    ncol = nsim
  )
  colnames(samples) <- paste0("sim", seq_len(nsim))

  class_pred <- NULL
  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0 || threshold >= 1) {
      stop("`threshold` must be a numeric value strictly between 0 and 1.")
    }
    class_pred <- ifelse(samples > threshold, 1L, 0L)
  }

  res <- list(
    samples = samples,
    mean = as.numeric(pred$mean),
    class = class_pred,
    X = object$X,
    Xnew = Xnew,
    threshold = threshold
  )
  class(res) <- "simulate_TwinBKP"
  res
}
