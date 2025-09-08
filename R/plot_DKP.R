#' @rdname plot
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
#'
#' # Plot results
#' plot(model1)
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
#'
#' # Plot results
#' plot(model2, n_grid = 50)
#'
#' @export
#' @method plot DKP

plot.DKP <- function(x, only_mean = FALSE, n_grid = 80, dims = NULL, ...){
  # ---------------- Argument Checking ----------------
  if (!is.logical(only_mean) || length(only_mean) != 1) {
    stop("`only_mean` must be a single logical value (TRUE or FALSE).")
  }

  if (!is.numeric(n_grid) || length(n_grid) != 1 || n_grid <= 0) {
    stop("'n_grid' must be a positive integer.")
  }
  n_grid <- as.integer(n_grid)

  # Extract necessary components from the DKP model object.
  X <- x$X # Covariate matrix.
  Y <- x$Y # Number of successes.
  Xbounds <- x$Xbounds

  d <- ncol(X)    # Dimensionality.
  q <- ncol(Y)    # Dimensionality.

  # Handle dims argument
  if (is.null(dims)) {
    if (d > 2) {
      stop("X has more than 2 dimensions. Please specify `dims` for plotting.")
    }
    dims <- seq_len(d)
  } else {
    if (!is.numeric(dims) || any(dims != as.integer(dims))) {
      stop("`dims` must be an integer vector.")
    }
    if (length(dims) < 1 || length(dims) > 2) {
      stop("`dims` must have length 1 or 2.")
    }
    if (any(dims < 1 | dims > d)) {
      stop(sprintf("`dims` must be within the range [1, %d].", d))
    }
    if (any(duplicated(dims))) {
      stop("`dims` cannot contain duplicate indices.")
    }
  }

  # Subset data to selected dimensions
  X_sub <- X[, dims, drop = FALSE]

  # old_par <- par(ask = TRUE)

  if (length(dims) == 1){
    #----- Plotting for 1-dimensional covariate data (d == 1) -----#

    # Generate new X values for smooth prediction
    Xnew <- matrix(seq(Xbounds[dims, 1], Xbounds[dims, 2], length.out = 10 * n_grid), ncol = 1)

    # Get the prediction for the new X values.
    Xnew_full <- lhs(nrow(Xnew), Xbounds)
    Xnew_full[, dims] <- Xnew
    prediction <- predict.DKP(x, Xnew_full, ...)

    # Determine whether it is a classification problem
    is_classification <- !is.null(prediction$class)

    old_par <- par(mfrow = c(2, 2))
    # on.exit(par(old_par))  # Restore par on exit

    # --- First panel: all mean curves together ---
    if(is_classification){
      cols <- rainbow(q)
      plot(NA, xlim = Xbounds, ylim = c(-0.1, 1.1),
           xlab = ifelse(d > 1, paste0("x", dims), "x"),
           ylab = "Probability",
           main = "Estimated Mean Curves (All Classes)")
      for (j in 1:q) {
        lines(Xnew, prediction$mean[, j], col = cols[j], lwd = 2)
      }
      for (i in 1:nrow(X)) {
        class_idx <- which.max(Y[i, ])
        points(X_sub[i], -0.05, col = cols[class_idx], pch = 20)
      }
      legend("top", legend = paste("Class", 1:q), col = cols, lty = 1, lwd = 2,
             horiz = TRUE, bty = "n")
    }

    # --- Remaining panels: each class with CI + obs ---
    for (j in 1:q) {
      mean_j  <- prediction$mean[, j]
      lower_j <- prediction$lower[, j]
      upper_j <- prediction$upper[, j]

      # Start plot for class j
      if (is_classification) {
        ylim = c(0, 1)
      }else{
        ylim = c(min(lower_j) * 0.9, min(1, max(upper_j) * 1.1))
      }
      plot(Xnew, mean_j,
           type = "l", col = "blue", lwd = 2,
           xlab = ifelse(d > 1, paste0("x", dims), "x"),
           ylab = "Probability",
           main = paste0("Estimated Probability (Class ", j, ")"),
           xlim = Xbounds,
           ylim = ylim)

      # Shaded CI
      polygon(c(Xnew, rev(Xnew)),
              c(lower_j, rev(upper_j)),
              col = "lightgrey", border = NA)
      lines(Xnew, mean_j, col = "blue", lwd = 2)

      # If class label is known, show binary observed indicator (1 if this class, 0 otherwise)
      if (is_classification) {
        obs_j <- as.integer(apply(Y, 1, which.max) == j)
        points(X_sub, obs_j, pch = 20, col = "red")
      } else {
        # Proportions from multinomial
        points(X_sub, Y[, j] / rowSums(Y), pch = 20, col = "red")

        # Legend
        if(j == 1) {
          legend("topleft",
                 legend = c("Estimated Probability",
                            paste0(prediction$CI_level * 100, "% CI"),
                            "Observed"),
                 col = c("blue", "lightgrey", "red"),
                 lwd = c(2, 8, NA), pch = c(NA, NA, 20), lty = c(1, 1, NA),
                 bty = "n")
        }
      }
    }
  } else {
    #----- Plotting for 2-dimensional covariate data (d == 2) -----#
    # Generate 2D prediction grid
    x1 <- seq(Xbounds[dims[1], 1], Xbounds[dims[1], 2], length.out = n_grid)
    x2 <- seq(Xbounds[dims[2], 1], Xbounds[dims[2], 2], length.out = n_grid)
    grid <- expand.grid(x1 = x1, x2 = x2)

    # Get the prediction for the new X values.
    Xnew_full <- lhs(nrow(grid), Xbounds)
    Xnew_full[, dims] <- as.matrix(grid)
    prediction <- predict.DKP(x, Xnew_full, ...)

    # Determine whether it is a classification problem
    is_classification <- !is.null(prediction$class)

    if(is_classification){
      df <- data.frame(x1 = grid$x1, x2 = grid$x2,
                       class = factor(prediction$class),
                       max_prob = apply(prediction$mean, 1, max))

      p1 <- my_2D_plot_fun_class("class", "Predicted Classes", df, X_sub, Y, dims= dims)
      p2 <- my_2D_plot_fun_class("max_prob", "Maximum Predicted Probability", df, X_sub, Y, classification = FALSE, dims= dims)
      grid.arrange(p1, p2, ncol = 2)
    }else{
      for (j in 1:q) {
        df <- data.frame(x1 = grid$x1, x2 = grid$x2,
                         Mean = prediction$mean[, j],
                         Upper = prediction$upper[, j],
                         Lower = prediction$lower[, j],
                         Variance = prediction$variance[, j])

        if (only_mean) {
          # Only plot the predicted mean graphs
          p1 <- my_2D_plot_fun("Mean", "Predictive Mean", df, dims= dims)
          print(p1)
        } else {
          # Create 4 plots
          p1 <- my_2D_plot_fun("Mean", "Predictive Mean", df, dims= dims)
          p2 <- my_2D_plot_fun("Upper", paste0(prediction$CI_level * 100, "% CI Upper"), df, dims= dims)
          p3 <- my_2D_plot_fun("Variance", "Predictive Variance", df, dims= dims)
          p4 <- my_2D_plot_fun("Lower", paste0(prediction$CI_level * 100, "% CI Lower"), df, dims= dims)
          # Arrange into 2×2 layout
          grid.arrange(p1, p2, p3, p4, ncol = 2,
                       top = textGrob(paste0("Estimated Probability (Class ", j, ")"),
                                      gp = gpar(fontface = "bold", fontsize = 16)))
        }
      }
    }
  }

  # par(old_par)
}
