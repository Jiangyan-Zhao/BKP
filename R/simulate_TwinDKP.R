#' @rdname simulate
#' @keywords TwinDKP
#' @export
#' @method simulate TwinDKP
simulate.TwinDKP <- function(object, nsim = 1, seed = NULL, Xnew = NULL, ...) {
  if (!is.numeric(nsim) || length(nsim) != 1 || nsim <= 0 || nsim != as.integer(nsim)) stop("`nsim` must be a positive integer.")
  if (!is.null(seed)) set.seed(seed)
  alpha_n <- if (is.null(Xnew)) object$alpha_n else predict(object, Xnew = Xnew, type = "probability", ...)$alpha_n
  n_new <- nrow(alpha_n); q <- ncol(alpha_n); samples <- array(0, c(n_new, q, as.integer(nsim)))
  for (i in seq_len(n_new)) samples[i,,] <- t(rdirichlet(as.integer(nsim), alpha_n[i, ]))
  dimnames(samples) <- list(paste0("x", seq_len(n_new)), paste0("Class", seq_len(q)), paste0("sim", seq_len(nsim)))
  out <- list(samples = samples, mean = alpha_n / rowSums(alpha_n), class = apply(samples, 3, max.col), X = object$X, Xnew = Xnew, ess = "none")
  class(out) <- "simulate_TwinDKP"; out
}
