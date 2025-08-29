#' @name print
#'
#' @title Print Summary of a Fitted BKP or DKP Model
#'
#' @description Displays a concise summary of a fitted BKP or DKP model. The
#'   output includes key characteristics such as sample size, input
#'   dimensionality, kernel type, loss function, optimized kernel
#'   hyperparameters, and minimum loss. Posterior predictive summaries (means
#'   and variances) are also provided.
#'
#' @param x An object of class \code{"BKP"} (from \code{\link{fit.BKP}}) or
#'   \code{"DKP"} (from \code{\link{fit.DKP}}).
#' @param ... Additional arguments passed to the generic \code{print} method
#'   (currently not used).
#'
#' @return Invisibly returns the input object (of class \code{"BKP"} or
#'   \code{"DKP"}). The function is called for its side effect of printing a
#'   summary to the console.
#'
#' @seealso \code{\link{fit.BKP}}, \code{\link{fit.DKP}},
#'   \code{\link{summary.BKP}}, \code{\link{summary.DKP}}.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv.
#'   https://doi.org/10.48550/arXiv.2508.10447.
#'
#' @keywords BKP
#'
#' @examples
#' # ============================================================== #
#' # ========================= BKP Examples ======================= #
#' # ============================================================== #
#'
#' #-------------------------- 1D Example ---------------------------
#' set.seed(123)
#'
#' # Define true success probability function
#' true_pi_fun <- function(x) {
#'   (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' }
#'
#' n <- 30
#' Xbounds <- matrix(c(-2,2), nrow=1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit BKP model
#' model1 <- fit.BKP(X, y, m, Xbounds=Xbounds)
#' print(model1)
#'
#'
#' #-------------------------- 2D Example ---------------------------
#' set.seed(123)
#'
#' # Define 2D latent function and probability transformation
#' true_pi_fun <- function(X) {
#'   if(is.null(nrow(X))) X <- matrix(X, nrow=1)
#'   m <- 8.6928
#'   s <- 2.4269
#'   x1 <- 4*X[,1]- 2
#'   x2 <- 4*X[,2]- 2
#'   a <- 1 + (x1 + x2 + 1)^2 *
#'     (19- 14*x1 + 3*x1^2- 14*x2 + 6*x1*x2 + 3*x2^2)
#'   b <- 30 + (2*x1- 3*x2)^2 *
#'     (18- 32*x1 + 12*x1^2 + 48*x2- 36*x1*x2 + 27*x2^2)
#'   f <- log(a*b)
#'   f <- (f- m)/s
#'   return(pnorm(f))  # Transform to probability
#' }
#'
#' n <- 100
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit BKP model
#' model2 <- fit.BKP(X, y, m, Xbounds=Xbounds)
#' print(model2)
#'
#' @export
#' @method print summary.BKP

print.summary.BKP <- function(x, ...) {
  cat("--------------------------------------------------\n")
  cat("        Summary of Beta Kernel Process (BKP)      \n")
  cat("--------------------------------------------------\n")
  cat(sprintf("Number of observations (n):  %d\n", x$n_obs))
  cat(sprintf("Input dimensionality (d):    %d\n", x$input_dim))
  cat(sprintf("Kernel type:                 %s\n", x$kernel))
  cat(sprintf("Optimized kernel parameters: %s\n",
              paste(sprintf("%.4f", x$theta_opt), collapse = ", ")))
  if (!is.na(x$loss_min)) {
    cat(sprintf("Minimum achieved loss:       %.5f\n", x$loss_min))
  }
  cat(sprintf("Loss function:               %s\n", x$loss))
  cat(sprintf("Prior type:                  %s\n", x$prior))
  if (x$prior == "fixed" || x$prior == "adaptive") {
    cat(sprintf("r0: %.3f\n", x$r0))
  }
  if (x$prior == "fixed") {
    cat(sprintf("p0: %.3f\n", x$p0))
  }

  # Posterior predictive summary
  cat("Posterior predictive summary (training points):\n")
  print(posterior_summary(x$post_mean, x$post_var))
  cat("--------------------------------------------------\n")

  invisible(x)
}


#' @rdname print
#'
#' @keywords BKP
#'
#' @export

print.BKP <- function(x, ...) {
  s <- summary.BKP(x)
  print(s)  # call print.summary.BKP
  invisible(x)
}


#' @rdname print
#'
#' @keywords BKP
#'
#' @export

print.predict.BKP <- function(x, ...) {
  n <- length(x$mean)

  # Determine prediction input
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
    X_preview <- data.frame(x1 = round(X_preview, 4))
    names(X_preview) <- "x"
  } else if (d == 2) {
    X_preview <- as.data.frame(round(X_preview, 4))
    names(X_preview) <- c("x1", "x2")
  } else {
    # Only display first and last columns with ... in between
    X_preview_vals <- round(X_preview[, c(1, d)], 3)
    X_preview <- as.data.frame(X_preview_vals)
    names(X_preview) <- c("x1", paste0("x", d))
    # Add a ... column
    X_preview$... <- rep("...", nrow(X_preview))
    # Reorder columns: x1, ..., xd
    X_preview <- X_preview[, c("x1", "...", paste0("x", d))]
  }

  # Construct results table
  pred_summary <- data.frame(
    mean     = round(head(x$mean, k), 4),
    variance = round(head(x$variance, k), 4),
    lower    = round(head(x$lower, k), 4),
    upper    = round(head(x$upper, k), 4)
  )

  # Update CI column names
  ci_low  <- round((1 - x$CI_level)/2 * 100, 2)
  ci_high <- round((1 + x$CI_level)/2 * 100, 2)
  names(pred_summary)[3:4] <- paste0(c(ci_low, ci_high), "% quantile")

  # Include predicted class if available
  if (!is.null(x$class)) {
    pred_summary$class <- head(x$class, k)
  }

  # Combine X preview and prediction
  res <- cbind(X_preview, pred_summary)

  print(res, row.names = FALSE)

  if (n > k) cat(" ...\n")

  invisible(x)
}
