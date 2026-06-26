#' @rdname print
#' @keywords TwinDKP
#' @export
#' @method print TwinDKP
print.TwinDKP <- function(x, ...) {
  cat("\n   Twin Dirichlet Kernel Process (TwinDKP) Model\n\n")
  cat(sprintf("Number of observations (n): %d\n", nrow(x$X)))
  cat(sprintf("Input dimension (d):        %d\n", ncol(x$X)))
  cat(sprintf("Number of classes (q):      %d\n", ncol(x$Y)))
  cat(sprintf("Global kernel:              %s\n", x$global_kernel)); cat(sprintf("Local kernel:               %s\n", x$local_kernel))
  cat(sprintf("Isotropic:                  %s\n", ifelse(x$isotropic, "TRUE", "FALSE")))
  cat(sprintf("theta_g:                    %s\n", paste(round(x$theta_g, 4), collapse = ", "))); cat(sprintf("theta_l:                    %.4f\n", x$theta_l))
  cat(sprintf("Loss function:              %s\n", x$loss)); cat(sprintf("Loss minimum:               %.5f\n", x$loss_min))
  cat(sprintf("Prior:                      %s\n", x$prior)); cat(sprintf("Global subset size (g):     %d (target %d)\n", length(x$global_indices), x$control$g_target))
  cat(sprintf("Local neighbours (l):       %d\n", x$control$l)); cat(sprintf("Twins runs:                 %d\n", x$control$twins)); invisible(x)
}
#' @rdname print
#' @keywords TwinDKP
#' @export
#' @method print predict_TwinDKP
print.predict_TwinDKP <- function(x, ...) { cat(if (is.null(x$Xnew)) "\nTwinDKP prediction on training data\n" else "\nTwinDKP prediction on new data\n"); print(head(x$mean)); invisible(x) }
#' @rdname print
#' @keywords TwinDKP
#' @export
#' @method print summary_TwinDKP
print.summary_TwinDKP <- function(x, ...) { cat("\n   Twin Dirichlet Kernel Process (TwinDKP) Model\n\n"); cat(sprintf("Number of observations (n): %d\n", x$n_obs)); cat(sprintf("Input dimension (d):        %d\n", x$input_dim)); cat(sprintf("Number of classes (q):      %d\n", x$n_classes)); cat(sprintf("Global subset size (g):     %d\n", x$global_size)); cat(sprintf("Local neighbours (l):       %d\n", x$local_size)); invisible(x) }
#' @rdname print
#' @keywords TwinDKP
#' @export
#' @method print simulate_TwinDKP
print.simulate_TwinDKP <- function(x, ...) { cat("\nTwinDKP posterior simulation\n"); print(dim(x$samples)); invisible(x) }
