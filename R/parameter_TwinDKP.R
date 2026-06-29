#' @rdname parameter
#' @keywords TwinDKP
#' @export
#' @method parameter TwinDKP
parameter.TwinDKP <- function(object, ...) {
  list(
    theta = object$theta_opt,
    theta_g = object$theta_g,
    theta_l = object$theta_l,
    alpha_n = object$alpha_n,
    global_kernel = object$global_kernel,
    local_kernel = object$local_kernel,
    global_indices = object$global_indices,
    control = object$control
  )
}
