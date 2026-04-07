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
  if (x$prior == "fixed" || x$prior == "adaptive") {
    cat(sprintf("r0: %.3f\n", x$r0))
  }
  if (x$prior == "fixed") {
    cat("p0: ", paste(round(x$p0, 3), collapse = ", "), "\n")
  }
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
  if (x$prior == "fixed" || x$prior == "adaptive") {
    cat(sprintf("    r0: %.3f\n", x$r0))
  }
  if (x$prior == "fixed") {
    cat("    p0: ", paste(round(x$p0, 3), collapse = ", "), "\n")
  }

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

  X_preview <- head(X_disp, k)
  if (d == 1) {
    X_preview <- data.frame(x1 = round(X_preview, 4))
    names(X_preview) <- "x"
  } else if (d == 2) {
    X_preview <- as.data.frame(round(X_preview, 4))
    names(X_preview) <- c("x1", "x2")
  } else {
    X_preview_vals <- round(X_preview[, c(1, d)], 3)
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

    res <- cbind(X_preview, pred_summary)
    print(res, row.names = FALSE)
    if (n > k) cat(" ...\n")
  }
  if (ncol(x$mean) > n_class) cat("\n ...\n")

  if (!is.null(x$class)) {
    cat("\nOverall predicted class (MAP):\n")
    print(head(x$class, k))
  }

  invisible(x)
}


#' @rdname print
#' @keywords DKP
#' @export
#' @method print simulate_TwinDKP
print.simulate_TwinDKP <- function(x, ...) {
  n <- dim(x$samples)[1]  # number of points
  q <- dim(x$samples)[2]  # number of classes
  nsim <- dim(x$samples)[3]  # number of simulations

  if (is.null(x$Xnew)) {
    cat("Simulation results on training data (X).\n")
    cat("Total number of training points:", n, "\n")
    X_disp <- x$X
  } else {
    cat("Simulation results on new data (Xnew).\n")
    cat("Total number of simulation points:", n, "\n")
    X_disp <- x$Xnew
  }
  cat("Number of posterior draws (nsim):", nsim, "\n")

  k <- min(6, n)
  d <- ncol(X_disp)

  if (n > k) {
    if (is.null(x$Xnew)) {
      cat("\nPreview of simulations for training data (first", k, "of", n, "points):\n")
    } else {
      cat("\nPreview of simulations for new data (first", k, "of", n, "points):\n")
    }
  } else {
    if (is.null(x$Xnew)) {
      cat("\nSimulations for all training data points:\n")
    } else {
      cat("\nSimulations for all new data points:\n")
    }
  }

  X_preview <- head(X_disp, k)
  if (d == 1) {
    X_preview <- data.frame(x1 = round(X_preview, 4))
    names(X_preview) <- "x"
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

  cat("\n--- Posterior Probability Simulations ---\n")

  show_sim <- min(3, nsim)
  for (s in seq_len(show_sim)) {
    cat("\nSimulation", s, ":\n")
    prob_mat <- as.matrix(x$samples[seq_len(k), , s])
    if (q <= 3) {
      prob_df <- round(prob_mat[, 1:q, drop = FALSE], 4)
      colnames(prob_df) <- paste0("Class", seq_len(q))
    } else {
      prob_df <- cbind(
        round(prob_mat[, 1:2, drop = FALSE], 4),
        "..." = rep("...", k),
        round(prob_mat[, q, drop = FALSE], 4)
      )
      colnames(prob_df)[c(1, 2, ncol(prob_df))] <- c("Class1", "Class2", paste0("Class", q))
    }

    print(cbind(X_preview, prob_df), row.names = FALSE)
    if (n > k) cat(" ...\n")
  }

  if (nsim > show_sim) cat("\nNote: only first 3 simulations are displayed out of", nsim, "simulations.\n")

  if (!is.null(x$class)) {
    class_preview <- head(x$class, k)
    if (nsim <= 3) {
      class_preview <- as.data.frame(class_preview)
      colnames(class_preview) <- paste0("sim", seq_len(nsim))
    } else {
      class_preview <- cbind(
        class_preview[, 1:2, drop = FALSE],
        "..." = rep("...", k),
        class_preview[, nsim, drop = FALSE]
      )
      colnames(class_preview)[c(1, 2, ncol(class_preview))] <- c("sim1", "sim2", paste0("sim", nsim))
    }
    cat("\n--- Classifications ---\n")
    print(cbind(X_preview, class_preview), row.names = FALSE)
    if (n > k) cat(" ...\n")
  }

  invisible(x)
}
