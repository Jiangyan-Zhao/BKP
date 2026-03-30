# Internal helper: compute posterior parameters on global support points
#' @keywords internal
.twinbkp_global_posterior <- function(object) {
  Xnorm_global <- object$Xnorm_global
  y_global     <- object$y_global
  m_global     <- object$m_global

  K_global <- kernel_matrix(
    X = Xnorm_global,
    theta = object$theta_global,
    kernel = object$kernel,
    isotropic = object$isotropic
  )

  prior_par <- get_prior(
    prior = object$prior,
    model = "BKP",
    r0 = object$r0,
    p0 = object$p0,
    y = y_global,
    m = m_global,
    K = K_global
  )

  post <- bkp_posterior_update_rcpp(
    K = K_global,
    y = as.numeric(y_global),
    m = as.numeric(m_global),
    alpha0 = as.numeric(prior_par$alpha0),
    beta0 = as.numeric(prior_par$beta0)
  )

  alpha_n <- as.numeric(post$alpha_n)
  beta_n  <- as.numeric(post$beta_n)
  mean_n  <- alpha_n / pmax(alpha_n + beta_n, 1e-10)
  var_n   <- mean_n * (1 - mean_n) / (alpha_n + beta_n + 1)

  list(
    alpha_n = alpha_n,
    beta_n = beta_n,
    mean = mean_n,
    variance = var_n
  )
}
