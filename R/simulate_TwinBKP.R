#' @rdname simulate
#'
#' @keywords BKP TwinBKP
#'
#' @export
#' @method simulate TwinBKP
simulate.TwinBKP <- function(object, nsim = 1, seed = NULL,
                             Xnew = NULL, threshold = NULL, ...) {
  if (!is.numeric(nsim) || length(nsim) != 1 || nsim <= 0 || nsim != as.integer(nsim)) {
    stop("`nsim` must be a positive integer.")
  }
  nsim <- as.integer(nsim)

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1 || seed != as.integer(seed))) {
    stop("`seed` must be a single integer or NULL.")
  }

  d <- ncol(object$X)
  if (!is.null(Xnew)) {
    if (is.null(dim(Xnew))) {
      if (d == 1L) {
        Xnew <- matrix(Xnew, ncol = 1L)
      } else {
        Xnew <- matrix(Xnew, nrow = 1L)
      }
    } else {
      Xnew <- as.matrix(Xnew)
    }
    if (!is.numeric(Xnew)) {
      stop("'Xnew' must be numeric.")
    }
    if (ncol(Xnew) != d) {
      stop("The number of columns in 'Xnew' must match the original input dimension.")
    }
  }

  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0 || threshold >= 1) {
      stop("`threshold` must be a numeric value strictly between 0 and 1 (e.g., 0.5).")
    }
  }

  if (!is.null(seed)) set.seed(seed)

  if (!is.null(Xnew)) {
    prediction <- predict.TwinBKP(object, Xnew = Xnew, type = "probability", ...)
    alpha_n <- prediction$alpha_n
    beta_n <- prediction$beta_n
  } else {
    alpha_n <- object$alpha_n
    beta_n <- object$beta_n
  }

  n_new <- if (!is.null(Xnew)) nrow(Xnew) else nrow(object$X)
  samples <- matrix(
    stats::rbeta(
      n_new * nsim,
      shape1 = rep(alpha_n, nsim),
      shape2 = rep(beta_n, nsim)
    ),
    nrow = n_new,
    ncol = nsim
  )
  colnames(samples) <- paste0("sim", 1:nsim)
  rownames(samples) <- paste0("x", 1:n_new)

  class_pred <- NULL
  if (!is.null(threshold)) {
    class_pred <- ifelse(samples > threshold, 1L, 0L)
  }

  simulation <- list(
    samples = samples,
    mean = alpha_n / (alpha_n + beta_n),
    class = class_pred,
    X = object$X,
    Xnew = Xnew,
    threshold = threshold
  )
  class(simulation) <- "simulate_TwinBKP"
  simulation
}
