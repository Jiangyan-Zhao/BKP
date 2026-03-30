#' @rdname print
#' @keywords DKP
#' @export
#' @method print TwinDKP
print.TwinDKP <- function(x, ...) {
  cat("\n       Twin Dirichlet Kernel Process (TwinDKP) Model   \n\n")
  cat(sprintf("Number of observations (n):  %d\n", nrow(x$X)))
  cat(sprintf("Input dimensionality (d):    %d\n", ncol(x$X)))
  cat(sprintf("Number of classes (q):       %d\n", ncol(x$Y)))
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
#' @keywords DKP
#' @export
#' @method print summary_TwinDKP
print.summary_TwinDKP <- function(x, ...) {
  cat("\n      Summary of Twin Dirichlet Kernel Process (TwinDKP) Model   \n\n")
  cat(sprintf("Number of observations (n):  %d\n", x$n_obs))
  cat(sprintf("Input dimensionality (d):    %d\n", x$input_dim))
  cat(sprintf("Number of classes (q):       %d\n", x$n_class))
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
  show_k <- min(3, x$n_class)
  for (j in seq_len(show_k)) {
    cat(sprintf("\nClass %d:\n", j))
    print(posterior_summary(x$post_mean_global[, j], x$post_var_global[, j]))
  }
  if (x$n_class > show_k) {
    cat("\n...\n")
    cat(sprintf("\nNote: Only the first %d classes are displayed out of %d classes.\n",
                show_k, x$n_class))
  }
  invisible(x)
}


#' @rdname print
#' @keywords DKP
#' @export
#' @method print predict_TwinDKP
print.predict_TwinDKP <- function(x, ...) {
  n <- nrow(x$mean)

  if (is.null(x$Xnew)) {
    cat("Prediction results on training data (X).\n")
    cat("Total number of training points:", n, "\n")
    X_disp <- x$X
  } else {
    cat("Prediction results on new data (Xnew).\n")
    cat("Total number of prediction points:", n, "\n")
    X_disp <- x$Xnew
  }

  d <- ncol(X_disp)
  k <- min(6, n)

  X_preview <- head(X_disp, k)
  if (d == 1) {
    X_preview <- data.frame(x = round(as.numeric(X_preview), 4))
  } else if (d == 2) {
    X_preview <- as.data.frame(round(X_preview, 4))
    names(X_preview) <- c("x1", "x2")
  } else {
    X_preview_vals <- round(X_preview[, c(1, d)], 4)
    X_preview <- as.data.frame(X_preview_vals)
    names(X_preview) <- c("x1", paste0("x", d))
    X_preview$... <- rep("...", nrow(X_preview))
    X_preview <- X_preview[, c("x1", "...", paste0("x", d))]
  }

  n_class <- min(3, ncol(x$mean))
  if (ncol(x$mean) > 3) {
    cat("\nNote: Only the first 3 classes are displayed out of", ncol(x$mean), "classes.\n")
  }

  ci_low  <- round((1 - x$CI_level)/2 * 100, 2)
  ci_high <- round((1 + x$CI_level)/2 * 100, 2)

  for (j in seq_len(n_class)) {
    cat("\nClass", j, "predictions:\n")
    pred_summary <- data.frame(
      Mean     = round(head(x$mean[, j], k), 4),
      Variance = round(head(x$variance[, j], k), 4),
      Lower    = round(head(x$lower[, j], k), 4),
      Upper    = round(head(x$upper[, j], k), 4)
    )
    names(pred_summary)[3:4] <- paste0(c(ci_low, ci_high), "% Quantile")
    print(cbind(X_preview, pred_summary), row.names = FALSE)
    if (n > k) cat(" ...\n")
  }
  if (ncol(x$mean) > n_class) cat("\n ...\n")

  if (!is.null(x$class)) {
    cat("\nOverall predicted class (MAP):\n")
    print(head(x$class, k))
  }

  invisible(x)
}
