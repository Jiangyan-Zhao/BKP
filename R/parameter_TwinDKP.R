#' @rdname parameter
#' @keywords TwinDKP
#' @export
#' @method parameter TwinDKP
parameter.TwinDKP <- function(object, ...) list(theta = object$theta_opt, theta_g = object$theta_g, theta_l = object$theta_l, prior = object$prior, r0 = object$r0, p0 = object$p0, global_kernel = object$global_kernel, local_kernel = object$local_kernel, alpha0 = object$alpha0, alpha_n = object$alpha_n, prob = fitted(object), global_indices = object$global_indices, control = object$control)
