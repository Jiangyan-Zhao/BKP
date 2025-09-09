#' @rdname print
#'
#' @keywords DKP
#'
#' @examples
#' # ============================================================== #
#' # ========================= DKP Examples ======================= #
#' # ============================================================== #
#'
#' #-------------------------- 1D Example ---------------------------
#' set.seed(123)
#'
#' # Define true class probability function (3-class)
#' true_pi_fun <- function(X) {
#'   p1 <- 1/(1+exp(-3*X))
#'   p2 <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
#'   return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
#' }
#'
#' n <- 30
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(150, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit DKP model
#' model1 <- fit_DKP(X, Y, Xbounds = Xbounds)
#' print(model1)                    # fitted object
#' print(summary(model1))           # summary
#' print(predict(model1))           # predictions
#' print(simulate(model1, nsim=3))  # posterior simulations
#'
#'
#' #-------------------------- 2D Example ---------------------------
#' set.seed(123)
#'
#' # Define latent function and transform to 3-class probabilities
#' true_pi_fun <- function(X) {
#'   if (is.null(nrow(X))) X <- matrix(X, nrow = 1)
#'   m <- 8.6928; s <- 2.4269
#'   x1 <- 4 * X[,1] - 2
#'   x2 <- 4 * X[,2] - 2
#'   a <- 1 + (x1 + x2 + 1)^2 *
#'     (19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2)
#'   b <- 30 + (2*x1 - 3*x2)^2 *
#'     (18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2)
#'   f <- (log(a*b)- m)/s
#'   p1 <- pnorm(f) # Transform to probability
#'   p2 <- sin(pi * X[,1]) * sin(pi * X[,2])
#'   return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
#' }
#'
#' n <- 100
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(150, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit DKP model
#' model2 <- fit_DKP(X, Y, Xbounds = Xbounds)
#' print(model2)                    # fitted object
#' print(summary(model2))           # summary
#' print(predict(model2))           # predictions
#' print(simulate(model2, nsim=3))  # posterior simulations
#'
#' @export
#' @method print DKP

print.DKP <- function(x, ...) {
  cat("\n       Dirichlet Kernel Process (DKP) Model   \n\n")
  cat(sprintf("Number of observations (n):  %d\n", nrow(x$X)))
  cat(sprintf("Input dimensionality (d):    %d\n", ncol(x$X)))
  cat(sprintf("Number of classes (q):       %d\n", ncol(x$Y)))
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
    cat("p0: ", paste(round(x$p0, 3), collapse = ", "), "\n")
  }

  invisible(x)
}


#' @rdname print
#'
#' @keywords DKP
#'
#' @export

print.summary_DKP <- function(x, ...) {
  cat("\n      Summary of Dirichlet Kernel Process (DKP) Model   \n\n")
  cat(sprintf("Number of observations (n):  %d\n", x$n_obs))
  cat(sprintf("Input dimensionality (d):    %d\n", x$input_dim))
  cat(sprintf("Number of classes (q):       %d\n", x$n_class))
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
    cat("    p0: ", paste(round(x$p0, 3), collapse = ", "), "\n")
  }

  # Posterior predictive summary
  n_class <- min(3, x$n_class)
  cat("\nPosterior predictive summary (training points):\n")

  for (j in 1:x$n_class) {
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
#'
#' @keywords DKP
#'
#' @export

print.predict_DKP <- function(x, ...) {
  n <- nrow(x$mean)

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

  # Only display first 3 classes or fewer
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

  # Add predicted class if available
  if (!is.null(x$class)) {
    cat("\nOverall predicted class (MAP):\n")
    print(head(x$class, k))
  }

  invisible(x)
}


#' @rdname print
#'
#' @keywords DKP
#'
#' @export

print.simulate_DKP <- function(x, ...) {
  n <- dim(x$samples)[1]  # number of points
  q <- dim(x$samples)[2]  # number of classes
  nsim <- dim(x$samples)[3]  # number of simulations

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

  # --- Posterior probabilities per simulation ---
  cat("\n--- Posterior Probability Simulations ---\n")

  nsim_to_show <- min(3, nsim)  # first 3 simulations
  for (i in 1:nsim_to_show) {
    cat("\nSimulation", i, ":\n")
    # Extract first 2 and last class
    prob_mat <- as.matrix(x$samples[1:k, , i])  # all classes for this sim
    if (q <= 3) {
      samples_preview <- round(prob_mat[, 1:q, drop = FALSE], 4)
      colnames(samples_preview) <- paste0("Class", 1:q)
    } else {
      samples_preview <- cbind(
        round(prob_mat[, 1:2, drop = FALSE], 4),
        "..." = rep("...", k),
        round(prob_mat[, q, drop = FALSE], 4)
      )
      colnames(samples_preview)[c(1,2,ncol(samples_preview))] <- c("Class1", "Class2", paste0("Class", q))
    }

    print(cbind(X_preview, samples_preview), row.names = FALSE)
    if (n > k) cat(" ...\n")
  }

  if (nsim > 3) cat("\nNote: only first 3 simulations are displayed out of", nsim, "simulations.\n")


  # If class predictions exist, include preview
  if (!is.null(x$class)) {
    class_preview <- head(x$class, k)
    if (nsim <= 3) {
      class_preview <- as.data.frame(class_preview)
      colnames(class_preview) <- paste0("sim", 1:nsim)
    } else {
      class_preview <- cbind(
        class_preview[, 1:2, drop = FALSE],
        "..." = rep("...", k),
        class_preview[, nsim, drop = FALSE]
      )
      colnames(class_preview)[c(1,2,ncol(class_preview))] <-
        c("sim1", "sim2", paste0("sim", nsim))
    }
    cat(paste0("\n--- Classifications ---\n"))
    print(cbind(X_preview, class_preview), row.names = FALSE)
    if (n > k) cat(" ...\n")
  }

  invisible(x)
}

