#' @rdname summary
#' @keywords TwinDKP
#' @export
#' @method summary TwinDKP
summary.TwinDKP <- function(object, ...) {
  alpha_n <- object$alpha_n
  row_sum <- rowSums(alpha_n)

  post_mean <- alpha_n / row_sum
  post_var <- (alpha_n * (row_sum - alpha_n)) / (row_sum^2 * (row_sum + 1))

  class_names <- paste0("class", seq_len(ncol(alpha_n)))
  colnames(post_mean) <- class_names
  colnames(post_var) <- class_names

  control <- object$control

  out <- list(
    n_obs = nrow(object$X),
    input_dim = ncol(object$X),
    n_class = ncol(object$Y),

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

  class(out) <- "summary_TwinDKP"
  return(out)
}
