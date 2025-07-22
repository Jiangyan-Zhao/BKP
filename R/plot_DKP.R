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
#'   p <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
#'   return(matrix(c(p/2, p/2, 1 - p), nrow = length(p)))
#' }
#'
#' n <- 30
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit DKP model
#' model1 <- fit.DKP(X, Y, Xbounds = Xbounds)
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
#'   f <- (log(a * b) - m) / s
#'   p <- pnorm(f)
#'   return(matrix(c(p/2, p/2, 1 - p), nrow = length(p)))
#' }
#'
#' n <- 100
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit DKP model
#' model2 <- fit.DKP(X, Y, Xbounds = Xbounds)
#'
#' # Plot results
#' plot(model2)
#'
#' @export
#' @method plot DKP

plot.DKP <- function(x, only_mean = FALSE, ...){
  if (!inherits(x, "DKP")) {
    stop("The input is not of class 'DKP'. Please provide a model fitted with 'fit.DKP()'.")
  }

  DKPmodel <- x

  # Extract necessary components from the DKP model object.
  X <- DKPmodel$X # Covariate matrix.
  Y <- DKPmodel$Y # Number of successes.
  Xbounds <- DKPmodel$Xbounds

  d <- ncol(X)    # Dimensionality.
  q <- ncol(Y)    # Dimensionality.

  if (d == 1){
    #----- Plotting for 1-dimensional covariate data (d == 1) -----#
    # Generate new X values for a smooth prediction curve.
    Xnew <- matrix(seq(Xbounds[1], Xbounds[2], length.out = 100), ncol = 1)

    # Get the prediction for the new X values.
    prediction <- predict.DKP(DKPmodel, Xnew, ...)

    for (j in 1:q) {
      mean_j <- prediction$mean[, j]
      lower_j <- prediction$lower[, j]
      upper_j <- prediction$upper[, j]
      # Initialize the plot with the estimated probability curve.
      plot(Xnew, mean_j,
           type = "l", col = "blue", lwd = 2,
           xlab = "x (Input Variable)", ylab = "Probability",
           main = paste0("Estimated Probability (class ", j, ")"),
           xlim = Xbounds,
           ylim = c(max(0, min(mean_j)-0.1),
                    min(1, max(mean_j)+0.2)))

      # Add a shaded Credible interval band using polygon.
      polygon(c(Xnew, rev(Xnew)),
              c(lower_j, rev(upper_j)),
              col = "lightgrey", border = NA)
      lines(Xnew, mean_j, col = "blue", lwd = 2)

      # Overlay observed proportions (y/m) as points.
      points(X, Y[,j] / rowSums(Y), pch = 20, col = "red")

      # Add a legend to explain plot elements.
      legend("topright",
             legend = c("Estimated Probability",
                        paste0((1 - prediction$CI_level)*100, "% Credible Interval"),
                        "Observed Proportions"),
             col = c("blue", "lightgrey", "red"), bty = "n",
             lwd = c(2, 8, NA), pch = c(NA, NA, 20), lty = c(1, 1, NA))
    }



  } else if (d == 2){
    #----- Plotting for 2-dimensional covariate data (d == 2) -----#
    # Generate 2D prediction grid
    x1 <- seq(Xbounds[1, 1], Xbounds[1, 2], length.out = 80)
    x2 <- seq(Xbounds[2, 1], Xbounds[2, 2], length.out = 80)
    grid <- expand.grid(x1 = x1, x2 = x2)
    prediction <- predict.DKP(DKPmodel, as.matrix(grid))

    for (j in 1:q) {
      df <- data.frame(x1 = grid$x1, x2 = grid$x2,
                       Mean = prediction$mean[, j],
                       Upper = prediction$upper[, j],
                       Lower = prediction$lower[, j],
                       Variance = prediction$variance[, j])

      if (only_mean) {
        # Only plot the predicted mean graphs
        my_2D_plot_fun("Mean", "Predictive Mean", df)
      } else {
        # Create 4 plots
        p1 <- my_2D_plot_fun("Mean", "Predictive Mean", df)
        p2 <- my_2D_plot_fun("Upper", paste0((1 - prediction$CI_level)*100, "% CI Upper"), df)
        p3 <- my_2D_plot_fun("Variance", "Predictive Variance", df)
        p4 <- my_2D_plot_fun("Lower", paste0((1 - prediction$CI_level)*100, "% CI Lower"), df)

        # Arrange into 2×2 layout
        grid.arrange(p1, p2, p3, p4, ncol = 2,
                     top = textGrob(paste0("Estimated Probability (class ", j, ")"),
                                    gp = gpar(fontface = "bold", fontsize = 16)))
      }
    }
  } else {
    # --- Error handling for higher dimensions ---
    stop("plot.DKP() only supports data where the dimensionality of X is 1 or 2.")
  }
}
