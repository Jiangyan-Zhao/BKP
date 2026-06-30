#' @rdname simulate
#' @keywords TwinDKP
#' @export
#' @method simulate TwinDKP
simulate.TwinDKP <- function(object, nsim = 1, seed = NULL,
                             Xnew = NULL, ...) {
  if (!is.numeric(nsim) || length(nsim) != 1L ||
      is.na(nsim) || !is.finite(nsim) ||
      nsim <= 0 || nsim != as.integer(nsim)) {
    stop("`nsim` must be a positive integer.")
  }
  nsim <- as.integer(nsim)

  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L ||
       is.na(seed) || !is.finite(seed) ||
       seed != as.integer(seed))) {
    stop("`seed` must be a single integer or NULL.")
  }

  if (!is.null(Xnew)) {
    prediction <- predict.TwinDKP(object, Xnew = Xnew, type = "probability", ...)
    alpha_n <- prediction$alpha_n
    Xnew <- prediction$Xnew
  } else {
    alpha_n <- object$alpha_n
  }

  n_new <- nrow(alpha_n)
  q <- ncol(alpha_n)
  samples <- array(0, c(n_new, q, as.integer(nsim)))
  for (i in seq_len(n_new)) {
    samples[i, , ] <- t(rdirichlet(as.integer(nsim), alpha_n[i, ]))
  }
  dimnames(samples) <- list(
    paste0("x", seq_len(n_new)),
    paste0("Class", seq_len(q)),
    paste0("sim", seq_len(nsim))
  )
  class_pred <- NULL
  if (all(rowSums(object$Y) == 1)) {
    class_pred <- matrix(NA_integer_, nrow = n_new, ncol = nsim)
    for (i in seq_len(nsim)) {
      class_pred[, i] <- max.col(samples[, , i])
    }
  }

  out <- list(
    samples = samples,
    mean = alpha_n / rowSums(alpha_n),
    class = class_pred,
    X = object$X,
    Xnew = Xnew
  )
  class(out) <- "simulate_TwinDKP"
  out
}
