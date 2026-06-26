#' @rdname predict
#' @keywords TwinDKP
#' @export
#' @method predict TwinDKP
predict.TwinDKP <- function(object, Xnew = NULL, CI_level = 0.95,
                            type = c("probability", "class", "count"),
                            Mnew = NULL, ...) {
  d <- ncol(object$X); type <- match.arg(type)
  if (!is.null(Xnew)) {
    if (is.null(dim(Xnew))) Xnew <- if (d == 1L) matrix(Xnew, ncol = 1L) else matrix(Xnew, nrow = 1L) else Xnew <- as.matrix(Xnew)
    if (!is.numeric(Xnew)) stop("'Xnew' must be numeric.")
    if (ncol(Xnew) != d) stop("The number of columns in 'Xnew' must match the original input dimension.")
    if (anyNA(Xnew) || any(!is.finite(Xnew))) stop("'Xnew' must contain only finite values with no NA, NaN, or Inf.")
  }
  n_pred <- if (is.null(Xnew)) nrow(object$X) else nrow(Xnew)
  if (!is.numeric(CI_level) || length(CI_level) != 1 || CI_level <= 0 || CI_level >= 1) stop("'CI_level' must be a single numeric value strictly between 0 and 1.")
  if (type == "count") {
    if (is.null(Mnew)) { if (is.null(Xnew)) Mnew <- rowSums(object$Y) else stop("When type = 'count' and Xnew is provided, 'Mnew' must also be provided.") }
    if (!is.numeric(Mnew) || !(length(Mnew) == 1L || length(Mnew) == n_pred) || anyNA(Mnew) || any(!is.finite(Mnew)) || any(Mnew <= 0) || any(Mnew != floor(Mnew))) stop("'Mnew' must contain positive integers and have length 1 or the same length as the number of prediction points.")
    Mnew <- as.integer(if (length(Mnew) == 1L) rep.int(Mnew, n_pred) else Mnew)
  }
  if (is.null(Xnew)) alpha_n <- object$alpha_n else {
    Xnew_norm <- sweep(sweep(Xnew, 2, object$Xbounds[,1], "-"), 2, object$Xbounds[,2] - object$Xbounds[,1], "/")
    local_indices <- twin_local_indices_rcpp(object$Xnorm, Xnew_norm, object$global_indices, object$control$l, if (is.null(object$control$leaf_size)) 8L else object$control$leaf_size)
    alpha_n <- .twindkp_compute_posterior(Xnew_norm, object$Xnorm, object$Y, object$global_indices, local_indices, object$theta_g, object$theta_l, object$global_kernel, object$local_kernel, object$isotropic, object$prior, object$r0, object$p0, FALSE)$alpha_n
  }
  alpha_n <- as.matrix(alpha_n); A <- rowSums(alpha_n); Amat <- matrix(A, nrow(alpha_n), ncol(alpha_n)); beta_n <- Amat - alpha_n
  prob_mean <- alpha_n / Amat; prob_var <- prob_mean * (1 - prob_mean) / (Amat + 1)
  class_names <- paste0("class", seq_len(ncol(alpha_n)))
  if (type == "class") return(structure(max.col(prob_mean), class = "predict_TwinDKP_class"))
  if (type == "probability") { pred_mean <- prob_mean; pred_var <- prob_var; pred_lower <- qbeta((1 - CI_level) / 2, alpha_n, beta_n); pred_upper <- qbeta((1 + CI_level) / 2, alpha_n, beta_n) } else {
    Mmat <- matrix(Mnew, nrow(alpha_n), ncol(alpha_n)); pred_mean <- Mmat * prob_mean; pred_var <- Mmat * (Amat + Mmat) * prob_var
    pred_lower <- matrix(qbetabinom_rcpp((1 - CI_level) / 2, as.vector(Mmat), as.vector(alpha_n), as.vector(beta_n)), nrow(alpha_n), ncol(alpha_n))
    pred_upper <- matrix(qbetabinom_rcpp((1 + CI_level) / 2, as.vector(Mmat), as.vector(alpha_n), as.vector(beta_n)), nrow(alpha_n), ncol(alpha_n))
  }
  colnames(pred_mean) <- colnames(pred_var) <- colnames(pred_lower) <- colnames(pred_upper) <- class_names
  out <- list(X = object$X, Xnew = Xnew, alpha_n = alpha_n, mean = pred_mean, variance = pred_var, lower = pred_lower, upper = pred_upper, CI_level = CI_level, type = type, class = max.col(prob_mean), ess = "none")
  if (type == "count") out$Mnew <- Mnew
  class(out) <- "predict_TwinDKP"; out
}
