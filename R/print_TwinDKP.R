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
  if (x$prior == "fixed" || x$prior == "adaptive") {
    cat(sprintf("r0:                         %.3f\n", x$r0))
  }
  if (x$prior == "fixed") {
    cat("p0:                         ", paste(round(x$p0, 3), collapse = ", "), "\n")
  }
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
  cat(sprintf("Number of classes (q):      %d\n", x$n_class))
  cat(sprintf("Global kernel:              %s\n", x$global_kernel))
  cat(sprintf("Local kernel:               %s\n", x$local_kernel))
  cat(sprintf("Isotropic:                  %s\n", ifelse(x$isotropic, "TRUE", "FALSE")))
  cat(sprintf("theta_g:                    %s\n", paste(round(x$theta_g, 4), collapse = ", ")))
  cat(sprintf("theta_l:                    %.4f\n", x$theta_l))
  cat(sprintf("Loss function:              %s\n", x$loss))
  cat(sprintf("Loss minimum:               %.5f\n", x$loss_min))
  cat(sprintf("Prior:                      %s\n", x$prior))

  if (x$prior == "fixed" || x$prior == "adaptive") {
    cat(sprintf("r0:                         %.3f\n", x$r0))
  }
  if (x$prior == "fixed") {
    cat("p0:                         ", paste(round(x$p0, 3), collapse = ", "), "\n")
  }

  cat(sprintf("Global subset size (g):     %d\n", x$global_size))
  cat(sprintf("Target global subset size:  %d\n", x$global_target))
  cat(sprintf("Local neighbours (l):       %d\n", x$local_size))
  cat(sprintf("Twinning runs:              %d\n", x$twins))

  n_class <- min(3, x$n_class)

  cat("\nPosterior predictive summary (training points):\n")

  for (j in seq_len(n_class)) {
    cat(sprintf("\nClass %d:\n", j))
    print(posterior_summary(x$post_mean[, j], x$post_var[, j]))
  }

  if (x$n_class > 3) {
    cat("\n ...\n")
    cat("\nNote: Only the first 3 classes are displayed out of", x$n_class, "classes.\n")
  }

  invisible(x)
}

#' @rdname print
#' @keywords TwinDKP
#' @export
#' @method print simulate_TwinDKP
print.simulate_TwinDKP <- function(x, ...) {
  n <- dim(x$samples)[1]
  q <- dim(x$samples)[2]
  nsim <- dim(x$samples)[3]

  # Determine simulation input
  if (is.null(x$Xnew)) {
    cat("TwinDKP simulation results on training data (X).\n")
    cat("Total number of training points:", n, "\n")
    X_disp <- x$X
  } else {
    cat("TwinDKP simulation results on new data (Xnew).\n")
    cat("Total number of simulation points:", n, "\n")
    X_disp <- x$Xnew
  }

  cat("Number of posterior draws (nsim):", nsim, "\n")

  d <- ncol(X_disp)
  k <- min(6, n)

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

  # Format X for display
  X_preview <- head(X_disp, k)
  if (d == 1) {
    X_preview <- data.frame(x = round(as.numeric(X_preview), 4))
  } else if (d == 2) {
    X_preview <- as.data.frame(round(X_preview, 4))
    names(X_preview) <- c("x1", "x2")
  } else {
    X_preview_vals <- round(X_preview[, c(1, d), drop = FALSE], 4)
    X_preview <- as.data.frame(X_preview_vals)
    names(X_preview) <- c("x1", paste0("x", d))
    X_preview$... <- rep("...", nrow(X_preview))
    X_preview <- X_preview[, c("x1", "...", paste0("x", d))]
  }

  # Posterior probabilities per simulation
  cat("\n--- TwinDKP Posterior Probability Simulations ---\n")

  nsim_to_show <- min(3, nsim)

  for (s in seq_len(nsim_to_show)) {
    cat("\nSimulation", s, ":\n")

    prob_mat <- matrix(
      x$samples[seq_len(k), , s],
      nrow = k,
      ncol = q
    )

    if (q <= 3) {
      samples_preview <- round(prob_mat[, seq_len(q), drop = FALSE], 4)
      colnames(samples_preview) <- paste0("Class", seq_len(q))
    } else {
      samples_preview <- cbind(
        round(prob_mat[, 1:2, drop = FALSE], 4),
        "..." = rep("...", k),
        round(prob_mat[, q, drop = FALSE], 4)
      )
      colnames(samples_preview)[c(1, 2, ncol(samples_preview))] <-
        c("Class1", "Class2", paste0("Class", q))
    }

    print(cbind(X_preview, samples_preview), row.names = FALSE)

    if (n > k) {
      cat(" ...\n")
    }
  }

  if (nsim > nsim_to_show) {
    cat("\nNote: only the first", nsim_to_show, "simulations are displayed out of",
        nsim, "simulations.\n")
  }

  # If class predictions exist, include preview
  if (!is.null(x$class)) {
    class_mat <- as.matrix(x$class)
    class_preview <- class_mat[seq_len(k), , drop = FALSE]

    if (nsim <= 3) {
      class_preview <- as.data.frame(class_preview)
      colnames(class_preview) <- paste0("sim", seq_len(nsim))
    } else {
      class_preview <- cbind(
        class_preview[, 1:2, drop = FALSE],
        "..." = rep("...", k),
        class_preview[, nsim, drop = FALSE]
      )
      colnames(class_preview)[c(1, 2, ncol(class_preview))] <-
        c("sim1", "sim2", paste0("sim", nsim))
    }

    cat("\n--- TwinDKP Classifications ---\n")
    print(cbind(X_preview, class_preview), row.names = FALSE)

    if (n > k) {
      cat(" ...\n")
    }
  }

  invisible(x)
}
