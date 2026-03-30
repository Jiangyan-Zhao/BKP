#' @rdname print
#' @keywords BKP
#' @export
#' @method print TwinBKP
print.TwinBKP <- function(x, ...) {
  cat("\n       Twin Beta Kernel Process (TwinBKP) Model    \n\n")
  cat(sprintf("Number of observations (n):  %d\n", nrow(x$X)))
  cat(sprintf("Input dimensionality (d):    %d\n", ncol(x$X)))
  cat(sprintf("Number of global points:     %d\n", length(x$global_idx)))
  cat(sprintf("g_nums:                      %d\n", x$g_nums))
  kernel_variant <- if (x$isotropic) "isotropic" else "anisotropic"
  cat(sprintf("Kernel type:                 (%s) %s\n", kernel_variant, x$kernel))
  cat(sprintf("Global theta:                %s\n",
              paste(sprintf("%.4f", x$theta_global), collapse = ", ")))
  cat(sprintf("Global loss:                 %.5f\n", x$loss_global))
  cat(sprintf("Loss function:               %s\n", x$loss))
  cat(sprintf("Prior type:                  %s\n", x$prior))
  invisible(x)
}


#' @rdname print
#' @keywords BKP
#' @export
#' @method print summary_TwinBKP
print.summary_TwinBKP <- function(x, ...) {
  cat("\n       Summary of Twin Beta Kernel Process (TwinBKP) Model   \n\n")
  cat(sprintf("Number of observations (n):  %d\n", x$n_obs))
  cat(sprintf("Input dimensionality (d):    %d\n", x$input_dim))
  cat(sprintf("Number of global points:     %d\n", x$n_global))
  cat(sprintf("g_nums:                      %d\n", x$g_nums))
  kernel_variant <- if (x$isotropic) "isotropic" else "anisotropic"
  cat(sprintf("Kernel type:                 (%s) %s\n", kernel_variant, x$kernel))
  cat(sprintf("Global theta:                %s\n",
              paste(sprintf("%.4f", x$theta_global), collapse = ", ")))
  cat(sprintf("Global loss:                 %.5f\n", x$loss_global))
  cat(sprintf("Loss function:               %s\n", x$loss))
  cat(sprintf("Prior type:                  %s\n", x$prior))

  cat("\nGlobal posterior summary (support points):\n")
  print(posterior_summary(x$post_mean_global, x$post_var_global))
  invisible(x)
}


#' @rdname print
#' @keywords BKP
#' @export
#' @method print simulate_TwinBKP
print.simulate_TwinBKP <- function(x, ...) {
  n <- nrow(x$samples)
  nsim <- ncol(x$samples)
  cat("Simulation results for TwinBKP.\n")
  cat("Total points:", n, "\n")
  cat("Number of posterior draws (nsim):", nsim, "\n")
  k <- min(6, n)
  out <- cbind(head(x$Xnew, k), round(head(x$samples, k), 4))
  print(out)
  if (n > k) cat("...\n")
  invisible(x)
}
