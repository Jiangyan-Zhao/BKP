#' @rdname summary
#' @keywords DKP
#'
#' @details
#' For \code{TwinDKP}, \code{summary()} focuses on global-stage tuning
#' information and posterior summaries computed on Twining support points.
#' Full prediction on arbitrary inputs should be obtained via
#' \code{predict.TwinDKP()}.
#'
#' @export
#' @method summary TwinDKP
summary.TwinDKP <- function(object, ...) {
  n_obs <- nrow(object$X)
  d <- ncol(object$X)
  q <- ncol(object$Y)
  if (is.null(object$mean_global) || is.null(object$var_global)) {
    stop("'mean_global'/'var_global' missing in object. Please refit using latest fit_TwinDKP().")
  }
  gp <- list(mean = object$mean_global, variance = object$var_global)

  res <- list(
    n_obs = n_obs,
    input_dim = d,
    n_class = q,
    n_global = length(object$global_idx),
    g_nums = object$g_nums,
    kernel = object$kernel,
    isotropic = object$isotropic,
    theta_global = object$theta_global,
    loss = object$loss,
    loss_global = object$loss_global,
    prior = object$prior,
    r0 = object$r0,
    p0 = object$p0,
    global_idx = object$global_idx,
    post_mean_global = gp$mean,
    post_var_global = gp$variance
  )

  class(res) <- "summary_TwinDKP"
  res
}
