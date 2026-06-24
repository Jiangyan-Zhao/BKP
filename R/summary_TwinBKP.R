#' @rdname summary
#'
#' @keywords BKP TwinBKP
#'
#' @export
#' @method summary TwinBKP
summary.TwinBKP <- function(object, ...) {

  n_obs <- nrow(object$X)
  d <- ncol(object$X)

  post_mean <- as.vector(object$alpha_n / (object$alpha_n + object$beta_n))
  post_var <- as.vector(
    post_mean * (1 - post_mean) / (object$alpha_n + object$beta_n + 1)
  )

  control <- object$control

  res <- list(
    n_obs = n_obs,
    input_dim = d,

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
    global_target = control$g_target,
    local_size = control$l,
    twins = control$twins,

    post_mean = post_mean,
    post_var = post_var
  )

  class(res) <- "summary_TwinBKP"
  res
}
