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

  if (ncol(X_disp) == 1) {
    X_preview <- data.frame(x = round(X_preview[,1], 4))
  } else {
    X_preview <- as.data.frame(round(X_preview, 4))
    names(X_preview) <- paste0("x", seq_len(ncol(X_preview)))
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
