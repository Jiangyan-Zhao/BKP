#' @rdname summary
#' @keywords TwinDKP
#' @export
#' @method summary TwinDKP
summary.TwinDKP <- function(object, ...) {
  prob <- fitted.TwinDKP(object)
  out <- list(
    n_obs = nrow(object$X),
    input_dim = ncol(object$X),
    n_classes = ncol(object$Y),
    kernel = object$kernel,
    global_kernel = object$global_kernel,
    local_kernel = object$local_kernel,
    isotropic = object$isotropic,
    theta_opt = object$theta_opt,
    theta_g = object$theta_g,
    theta_l = object$theta_l,
    loss = object$loss,
    loss_min = object$loss_min,
    prior = object$prior,
    r0 = object$r0,
    p0 = object$p0,
    global_size = length(object$global_indices),
    global_target = object$control$g_target,
    local_size = object$control$l,
    twins = object$control$twins,
    post_mean = prob)
  class(out) <- "summary_TwinDKP"
  out
}
