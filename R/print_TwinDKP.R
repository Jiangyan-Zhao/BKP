#' @rdname print
#' @keywords TwinDKP
#' @export
#' @method print TwinDKP
print.TwinDKP <- function(x, ...) {
  cat("\n   Twin Dirichlet Kernel Process (TwinDKP) Model\n\n")
  cat(sprintf("Number of observations (n): %d\n", nrow(x$X)))
  cat(sprintf("Input dimension (d):        %d\n", ncol(x$X)))
  cat(sprintf("Number of classes (q):      %d\n", ncol(x$Y)))
  cat(sprintf("Global kernel:              %s\n", x$global_kernel))
  cat(sprintf("Local kernel:               %s\n", x$local_kernel))
  cat(sprintf("Isotropic:                  %s\n", ifelse(x$isotropic, "TRUE", "FALSE")))
  cat(sprintf("theta_g:                    %s\n", paste(round(x$theta_g, 4), collapse = ", ")))
  cat(sprintf("theta_l:                    %.4f\n", x$theta_l))
  cat(sprintf("Loss function:              %s\n", x$loss))
  cat(sprintf("Loss minimum:               %.5f\n", x$loss_min))
  cat(sprintf("Prior:                      %s\n", x$prior))
  cat(sprintf(
    "Global subset size (g):     %d (target %d)\n",
    length(x$global_indices), x$control$g_target
  ))
  cat(sprintf("Local neighbours (l):       %d\n", x$control$l))
  cat(sprintf("Twins runs:                 %d\n", x$control$twins))
  invisible(x)
}

#' @rdname print
#' @keywords TwinDKP
#' @export
#' @method print predict_TwinDKP
print.predict_TwinDKP <- function(x, ...) {
  n <- nrow(x$mean)

  # Determine prediction input
  if (is.null(x$Xnew)) {
    cat("TwinDKP prediction results on training data (X).\n")
    cat("Total number of training points:", n, "\n")
    X_disp <- x$X
  } else {
    cat("TwinDKP prediction results on new data (Xnew).\n")
    cat("Total number of prediction points:", n, "\n")
    X_disp <- x$Xnew
  }

  d <- ncol(X_disp)

  # Determine how many rows to preview
  k <- min(6, n)
  if (n > k) {
    if (is.null(x$Xnew)) {
      cat("\nPreview of predictions for training data (first", k, "of", n, "points):\n")
    } else {
      cat("\nPreview of predictions for new data (first", k, "of", n, "points):\n")
    }
  } else {
    if (is.null(x$Xnew)) {
      cat("\nPredictions for all training data points:\n")
    } else {
      cat("\nPredictions for all new data points:\n")
    }
  }

  # Format X_disp for display
  X_preview <- head(X_disp, k)
  if (d == 1) {
    X_preview <- data.frame(x = round(as.numeric(X_preview), 4))
  } else if (d == 2) {
    X_preview <- as.data.frame(round(X_preview, 4))
    names(X_preview) <- c("x1", "x2")
  } else {
    X_preview_vals <- round(X_preview[, c(1, d), drop = FALSE], 3)
    X_preview <- as.data.frame(X_preview_vals)
    names(X_preview) <- c("x1", paste0("x", d))
    X_preview$... <- rep("...", nrow(X_preview))
    X_preview <- X_preview[, c("x1", "...", paste0("x", d))]
  }

  # Only display first 3 classes or fewer
  n_class <- min(3, ncol(x$mean))
  if (ncol(x$mean) > 3) {
    cat("\nNote: Only the first 3 classes are displayed out of", ncol(x$mean), "classes.\n")
  }

  ci_low <- round((1 - x$CI_level) / 2 * 100, 2)
  ci_high <- round((1 + x$CI_level) / 2 * 100, 2)

  for (j in seq_len(n_class)) {
    cat("\nClass", j, "predictions:\n")
    pred_summary <- data.frame(
      Mean = round(head(x$mean[, j], k), 4),
      Variance = round(head(x$variance[, j], k), 4),
      Lower = round(head(x$lower[, j], k), 4),
      Upper = round(head(x$upper[, j], k), 4)
    )
    names(pred_summary)[3:4] <- paste0(c(ci_low, ci_high), "% Quantile")

    res <- cbind(X_preview, pred_summary)
    print(res, row.names = FALSE)

    if (n > k) {
      cat(" ...\n")
    }
  }

  if (ncol(x$mean) > n_class) {
    cat("\n ...\n")
  }

  if (!is.null(x$class)) {
    cat("\nOverall predicted class (MAP):\n")
    print(head(x$class, k))
  }

  invisible(x)
}

#' @rdname print
#' @keywords TwinDKP
#' @export
#' @method print summary_TwinDKP
print.summary_TwinDKP <- function(x, ...) {
  cat("\n   Twin Dirichlet Kernel Process (TwinDKP) Model\n\n")
  cat(sprintf("Number of observations (n): %d\n", x$n_obs))
  cat(sprintf("Input dimension (d):        %d\n", x$input_dim))
  cat(sprintf("Number of classes (q):      %d\n", x$n_classes))
  cat(sprintf("Global subset size (g):     %d\n", x$global_size))
  cat(sprintf("Local neighbours (l):       %d\n", x$local_size))
  invisible(x)
}

#' @rdname print
#' @keywords TwinDKP
#' @export
#' @method print simulate_TwinDKP
print.simulate_TwinDKP <- function(x, ...) {
  cat("\nTwinDKP posterior simulation\n")
  print(dim(x$samples))
  invisible(x)
}
