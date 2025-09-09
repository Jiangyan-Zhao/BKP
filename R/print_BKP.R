#' @name print
#'
#' @title Print Methods for BKP and DKP Objects
#'
#' @description Provides formatted console output for fitted BKP/DKP model
#'   objects, their summaries, predictions, and simulations. The following
#'   specialized methods are supported:
#'   \itemize{
#'     \item \code{print.BKP}, \code{print.DKP} – display fitted model objects.
#'     \item \code{print.summary_BKP}, \code{print.summary_DKP} – display
#'       model summaries.
#'     \item \code{print.predict_BKP}, \code{print.predict_DKP} – display
#'       posterior predictive results.
#'     \item \code{print.simulate_BKP}, \code{print.simulate_DKP} – display
#'       posterior simulations.
#'   }
#'
#' @param x An object of class \code{"BKP"} or \code{"DKP"}, or a derived object
#'   such as \code{summary}, \code{predict}, or \code{simulate}.
#' @param ... Additional arguments passed to the generic \code{print} method
#'   (currently unused; included for S3 consistency).
#'
#' @return Invisibly returns the input object. Called for the side effect of
#'   printing human-readable summaries to the console.
#'
#' @seealso \code{\link{fit_BKP}}, \code{\link{fit_DKP}} for model fitting;
#'   \code{\link{summary.BKP}}, \code{\link{summary.DKP}} for model summaries;
#'   \code{\link{predict.BKP}}, \code{\link{predict.DKP}} for posterior
#'   prediction; \code{\link{simulate.BKP}}, \code{\link{simulate.DKP}} for
#'   posterior simulations.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv.
#'   https://doi.org/10.48550/arXiv.2508.10447.
#'
#' @keywords BKP DKP
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
#' model1 <- fit_BKP(X, y, m, Xbounds=Xbounds)
#' print(model1)                    # fitted object
#' print(summary(model1))           # summary
#' print(predict(model1))           # predictions
#' print(simulate(model1, nsim=3))  # posterior simulations
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
#' model2 <- fit_BKP(X, y, m, Xbounds=Xbounds)
#' print(model2)                    # fitted object
#' print(summary(model2))           # summary
#' print(predict(model2))           # predictions
#' print(simulate(model2, nsim=3))  # posterior simulations
#'
#' @export
#' @method print BKP

print.BKP <- function(x, ...) {
  cat("\n       Beta Kernel Process (BKP) Model    \n\n")
  cat(sprintf("Number of observations (n):  %d\n", nrow(x$X)))
  cat(sprintf("Input dimensionality (d):    %d\n", ncol(x$X)))
  cat(sprintf("Kernel type:                 %s\n", x$kernel))
  cat(sprintf("Optimized kernel parameters: %s\n",
              paste(sprintf("%.4f", x$theta_opt), collapse = ", ")))
  cat(sprintf("Minimum achieved loss:       %.5f\n", x$loss_min))
  cat(sprintf("Loss function:               %s\n", x$loss))
  cat(sprintf("Prior type:                  %s\n", x$prior))
  if (x$prior == "fixed" || x$prior == "adaptive") {
    cat(sprintf("r0: %.3f\n", x$r0))
  }
  if (x$prior == "fixed") {
    cat(sprintf("p0: %.3f\n", x$p0))
  }

  invisible(x)
}


#' @rdname print
#'
#' @keywords BKP
#'
#' @export

print.summary_BKP <- function(x, ...) {
  cat("\n       Beta Kernel Process (BKP) Model   \n\n")
  cat(sprintf("Number of observations (n):  %d\n", x$n_obs))
  cat(sprintf("Input dimensionality (d):    %d\n", x$input_dim))
  cat(sprintf("Kernel type:                 %s\n", x$kernel))
  cat(sprintf("Optimized kernel parameters: %s\n",
              paste(sprintf("%.4f", x$theta_opt), collapse = ", ")))
  cat(sprintf("Minimum achieved loss:       %.5f\n", x$loss_min))
  cat(sprintf("Loss function:               %s\n", x$loss))
  cat(sprintf("Prior type:                  %s\n", x$prior))
  if (x$prior == "fixed" || x$prior == "adaptive") {
    cat(sprintf("    r0: %.3f\n", x$r0))
  }
  if (x$prior == "fixed") {
    cat(sprintf("    p0: %.3f\n", x$p0))
  }

  # Posterior predictive summary
  cat("\nPosterior predictive summary (training points):\n")
  print(posterior_summary(x$post_mean, x$post_var))

  invisible(x)
}


#' @rdname print
#'
#' @keywords BKP
#'
#' @export

print.predict_BKP <- function(x, ...) {
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
    X_preview_vals <- round(X_preview[, c(1, d)], 4)
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


#' @rdname print
#'
#' @keywords BKP
#'
#' @export

print.simulate_BKP <- function(x, ...) {
  n <- length(x$mean)
  nsim <- ncol(x$samples)

  # Determine simulation input
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

  # Determine how many rows to preview
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
    X_preview <- data.frame(x1 = round(X_preview, 4))
    names(X_preview) <- "x"
  } else if (d == 2) {
    X_preview <- as.data.frame(round(X_preview, 4))
    names(X_preview) <- c("x1", "x2")
  } else {
    # Only display first and last columns with ... in between
    X_preview_vals <- round(X_preview[, c(1, d)], 4)
    X_preview <- as.data.frame(X_preview_vals)
    names(X_preview) <- c("x1", paste0("x", d))
    # Add a ... column
    X_preview$... <- rep("...", nrow(X_preview))
    # Reorder columns: x1, ..., xd
    X_preview <- X_preview[, c("x1", "...", paste0("x", d))]
  }

  # Add first few simulation columns
  samples_preview <- head(x$samples, k)
  if (nsim <= 3) {
    samples_preview <- round(samples_preview, 4)
    samples_preview <- as.data.frame(samples_preview)
  } else {
    # show first 2 and last sim with "..."
    samples_preview <- cbind(
      round(samples_preview[, 1:2, drop = FALSE], 4),
      "..." = rep("...", k),
      round(samples_preview[, nsim, drop = FALSE], 4)
    )
    colnames(samples_preview)[c(1,2,ncol(samples_preview))] <-
      c("sim1", "sim2", paste0("sim", nsim))
  }

  cat("\n--- Posterior Probability Simulations ---\n")
  print(cbind(X_preview, samples_preview), row.names = FALSE)
  if (n > k) cat(" ...\n")

  # If class predictions exist, include preview
  if (!is.null(x$class)) {
    class_preview <- head(x$class, k)
    if (nsim <= 3) {
      class_preview <- as.data.frame(class_preview)
    } else {
      class_preview <- cbind(
        class_preview[, 1:2, drop = FALSE],
        "..." = rep("...", k),
        class_preview[, nsim, drop = FALSE]
      )
      colnames(class_preview)[c(1,2,ncol(class_preview))] <-
        c("sim1", "sim2", paste0("sim", nsim))
    }
    cat(paste0("\n--- Binary Classifications (threshold = ", x$threshold, ") ---\n"))
    print(cbind(X_preview, class_preview), row.names = FALSE)
    if (n > k) cat(" ...\n")
  }

  invisible(x)
}
