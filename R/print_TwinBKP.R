#' @rdname print
#'
#' @keywords BKP TwinBKP
#'
#' @export
#' @method print TwinBKP
print.TwinBKP <- function(x, ...) {

  cat("\n   Twin Beta Kernel Process (TwinBKP) Model\n\n")

  cat(sprintf("Number of observations (n): %d\n", nrow(x$X)))
  cat(sprintf("Input dimension (d):        %d\n", ncol(x$X)))

  cat(sprintf("Global kernel:              %s\n", x$global_kernel))
  cat(sprintf("Local kernel:               %s\n", x$local_kernel))

  cat(sprintf("Isotropic:                  %s\n",
              ifelse(x$isotropic, "TRUE", "FALSE")))

  cat(sprintf("theta_g:                    %s\n",
              paste(round(x$theta_g, 4), collapse = ", ")))
  cat(sprintf("theta_l:                    %.4f\n", x$theta_l))

  cat(sprintf("Loss function:              %s\n", x$loss))
  cat(sprintf("Loss minimum:               %.5f\n", x$loss_min))

  cat(sprintf("Prior:                      %s\n", x$prior))

  cat(sprintf("Global subset size (g):     %d (target %d)\n",
              length(x$global_indices),
              x$control$g_target))

  cat(sprintf("Local neighbours (l):       %d\n", x$control$l))
  cat(sprintf("Twins runs:                 %d\n", x$control$twins))

  invisible(x)
}


#' @rdname print
#'
#' @keywords BKP TwinBKP
#'
#' @export
#' @method print predict_TwinBKP
print.predict_TwinBKP <- function(x, ...) {

  n <- length(x$mean)

  if (is.null(x$Xnew)) {
    cat("\nTwinBKP prediction on training data\n")
    X_disp <- x$X
  } else {
    cat("\nTwinBKP prediction on new data\n")
    X_disp <- x$Xnew
  }

  cat(sprintf("Number of points: %d\n\n", n))

  k <- min(6, n)
  X_preview <- head(X_disp, k)

  d <- ncol(X_disp)
  if (d == 1) {
    X_preview <- data.frame(x = round(X_preview[, 1], 4))
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

  pred_summary <- data.frame(
    mean = round(head(x$mean, k), 4),
    variance = round(head(x$variance, k), 4),
    lower = round(head(x$lower, k), 4),
    upper = round(head(x$upper, k), 4)
  )

  ci_low <- round((1 - x$CI_level)/2 * 100, 1)
  ci_high <- round((1 + x$CI_level)/2 * 100, 1)

  names(pred_summary)[3:4] <- c(
    paste0(ci_low, "%"),
    paste0(ci_high, "%")
  )

  if (!is.null(x$class)) {
    pred_summary$class <- head(x$class, k)
  }

  out <- cbind(X_preview, pred_summary)
  print(out, row.names = FALSE)

  invisible(x)
}

#' @rdname print
#'
#' @keywords BKP TwinBKP
#'
#' @export
#' @method print summary_TwinBKP
print.summary_TwinBKP <- function(x, ...) {
  cat("\n   Twin Beta Kernel Process (TwinBKP) Model\n\n")

  cat(sprintf("Number of observations (n): %d\n", x$n_obs))
  cat(sprintf("Input dimension (d):        %d\n", x$input_dim))
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
    cat(sprintf("p0:                         %.3f\n", x$p0))
  }
  cat(sprintf("Global subset size (g):     %d\n", x$global_size))
  cat(sprintf("Target global subset size:  %d\n", x$global_target))
  cat(sprintf("Local neighbours (l):       %d\n", x$local_size))
  cat(sprintf("Twinning runs:              %d\n", x$twins))

  cat("\nPosterior predictive summary (training points):\n")
  print(posterior_summary(x$post_mean, x$post_var))

  invisible(x)
}

#' @rdname print
#'
#' @keywords BKP TwinBKP
#'
#' @export
#' @method print simulate_TwinBKP
print.simulate_TwinBKP <- function(x, ...) {
  n <- length(x$mean)
  nsim <- ncol(x$samples)

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

  X_preview <- head(X_disp, k)
  if (d == 1) {
    X_preview <- data.frame(x = round(X_preview[, 1], 4))
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

  samples_preview <- head(x$samples, k)
  if (nsim <= 3) {
    samples_preview <- as.data.frame(round(samples_preview, 4))
  } else {
    samples_preview <- cbind(
      round(samples_preview[, 1:2, drop = FALSE], 4),
      "..." = rep("...", k),
      round(samples_preview[, nsim, drop = FALSE], 4)
    )
    colnames(samples_preview)[c(1, 2, ncol(samples_preview))] <-
      c("sim1", "sim2", paste0("sim", nsim))
  }

  cat("\n--- TwinBKP Posterior Probability Simulations ---\n")
  print(cbind(X_preview, samples_preview), row.names = FALSE)
  if (n > k) cat(" ...\n")

  if (!is.null(x$class)) {
    cat(paste0("\nBinary class labels available (threshold = ", x$threshold, ").\n"))
  }

  invisible(x)
}
