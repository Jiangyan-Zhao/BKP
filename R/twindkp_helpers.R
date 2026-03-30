# Internal helper: compute posterior parameters on global support points
#' @keywords internal
.twindkp_global_posterior <- function(object) {
  Xnorm_global <- object$Xnorm_global
  Y_global     <- object$Y_global

  K_global <- kernel_matrix(
    X = Xnorm_global,
    theta = object$theta_global,
    kernel = object$kernel,
    isotropic = object$isotropic
  )

  alpha0 <- get_prior(
    prior = object$prior,
    model = "DKP",
    r0 = object$r0,
    p0 = object$p0,
    Y = Y_global,
    K = K_global
  )

  alpha_n <- as.matrix(alpha0) + as.matrix(K_global %*% Y_global)
  row_sum <- rowSums(alpha_n)
  mean_n  <- alpha_n / pmax(row_sum, 1e-10)
  var_n   <- mean_n * (1 - mean_n) / (row_sum + 1)

  class_names <- paste0("class", seq_len(ncol(alpha_n)))
  colnames(alpha_n) <- class_names
  colnames(mean_n) <- class_names
  colnames(var_n) <- class_names

  list(
    alpha_n = alpha_n,
    mean = mean_n,
    variance = var_n
  )
}
