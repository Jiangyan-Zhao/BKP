#' @name plot
#'
#' @title Plot a Beta Kernel Process (BKP) Model
#'
#' @description
#' Visualizes the fitted Beta Kernel Process (BKP) model according to the dimensionality of the covariate matrix `X`.
#' For 1-dimensional input, the function plots the estimated probability curve with confidence intervals and observed data points.
#' For 2-dimensional input, it produces contour plots for the estimated mean, confidence bounds, and uncertainty.
#'
#' @param x An object of class \code{"BKP"}, typically returned by the \code{\link{fit.BKP}} function.
#' @param ... Additional arguments passed to generic plot functions (currently not used, included for S3 method consistency).
#'
#' @details
#' The plotting behavior depends on the dimensionality of the covariates \code{X} in the fitted BKP model:
#' \itemize{
#'   \item \strong{1D covariates:} Produces a line plot of the predicted probability across the input domain.
#'         A shaded region represents the confidence interval (default 95%), and the observed proportions \code{y/m}
#'         are added as points.
#'   \item \strong{2D covariates:} Generates a 2x2 grid of contour plots showing the predicted mean, upper bound,
#'         lower bound, and the width of the confidence interval (uncertainty).
#' }
#'
#' If the dimension of \code{X} exceeds 2, the function will stop with an informative error.
#'
#' @author Jiangyan Zhao, Kunhai Qing, Jin Xu
#'
#' @seealso
#' \code{\link{fit.BKP}} for fitting a BKP model.
#' \code{\link{predict.BKP}} for generating predictions from a BKP model.
#'
#' @keywords BKP
#'
#' @export
#' @method plot BKP



plot.BKP <- function(x, ...){
  if (!inherits(x, "BKP")) {
    stop("The input is not of class 'BKP'. Please provide a model fitted with 'fit.BKP()'.")
  }

  BKPmodel <- x

  # Extract necessary components from the BKP model object.
  X <- BKPmodel$X # Covariate matrix.
  y <- BKPmodel$y # Number of successes.
  m <- BKPmodel$m # Number of trials.
  Xbounds <- BKPmodel$Xbounds

  d <- ncol(X)    # Dimensionality.

  if (d == 1){
    #----- Plotting for 1-dimensional covariate data (d == 1) -----#
    # Generate new X values for a smooth prediction curve.
    Xnew <- matrix(seq(Xbounds[1], Xbounds[2], length.out = 1000), ncol = 1)

    # Get the prediction for the new X values.
    prediction <- predict(BKPmodel, Xnew, CI.size = 0.05)

    # Initialize the plot with the estimated probability curve.
    plot(prediction$X, prediction$mean,
         type = "l", col = "blue", lwd = 2,
         xlab = "x (Input Variable)", ylab = "Probability",
         main = "BKP Estimated Probability",
         xlim = Xbounds,
         ylim = c(max(0, min(prediction$mean)-0.1),
                  min(1, max(prediction$mean)+0.1)))

    # Add a shaded confidence interval band using polygon.
    polygon(c(prediction$X, rev(prediction$X)),
            c(prediction$lower, rev(prediction$upper)),
            col = "lightgrey", border = NA)
    lines(prediction$X, prediction$mean, col = "blue", lwd = 2)

    # Overlay observed proportions (y/m) as points.
    points(X, y / m, pch = 16, cex = 0.8, col = "red")

    # Add a legend to explain plot elements.
    legend("topright",
           legend = c("Estimated Probability", "95% Confidence Interval", "Observed Proportions"),
           col = c("blue", "lightgrey", "red"),
           lwd = c(2, NA, NA), pch = c(NA, 15, 16), lty = c(1, NA, NA), bty = "n")
  } else if (d == 2){
    #----- Plotting for 2-dimensional covariate data (d == 2) -----#
    # Create a grid of new X values for prediction over the 2D space.
    x1 <- seq(Xbounds[1,1], Xbounds[1,2], length.out = 100)
    x2 <- seq(Xbounds[2,1], Xbounds[2,2], length.out = 100)
    Xnew <- expand.grid(x1 = x1, x2 = x2)
    Xnew <- as.matrix(Xnew)

    # Predict probabilities and confidence intervals for the grid.
    prediction <- predict(BKPmodel, Xnew, CI.size = 0.05)

    # Reshape the 1D prediction results into 2D matrices for image plotting.
    pred_mean <- matrix(prediction$mean, nrow = length(x1), ncol = length(x2))
    pred_upper <- matrix(prediction$upper, nrow = length(x1), ncol = length(x2))
    pred_lower <- matrix(prediction$lower, nrow = length(x1), ncol = length(x2))
    pred_variance <- matrix(prediction$variance, nrow = length(x1), ncol = length(x2))

    # Set up a 2x2 plot layout for multiple contour plots.
    layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))

    # Plot 1: Estimated Probability Contour Map
    image(x1, x2, pred_mean,
          xlab = "X1", ylab = "X2",
          main = "Mean",
          col = hcl.colors(100, "viridis"))
    contour(x1, x2, pred_mean, add = TRUE, col = "black")
    # Overlay observed proportions (y/m) as points.

    # Plot 2: Upper Bound Contour Map
    image(x1, x2, pred_upper,
          xlab = "X1", ylab = "X2",
          main = "Upper Bound",
          col = hcl.colors(100, "viridis"))
    contour(x1, x2, pred_upper, add = TRUE, col = "black")

    # Plot 3: Lower Bound Contour Map
    image(x1, x2, pred_lower,
          xlab = "X1", ylab = "X2",
          main = "Lower Bound",
          col = hcl.colors(100, "viridis"))
    contour(x1, x2, pred_lower, add = TRUE, col = "black")

    # Plot 4: Variance Contour Map
    image(x1, x2, pred_variance,
          xlab = "X1", ylab = "X2",
          main = "Variance",
          col = hcl.colors(100, "plasma"))
    contour(x1, x2, pred_variance, add = TRUE, col = "black")

    par(mfrow = c(1, 1)) # Reset the plot layout to default (1 plot per page) after generating the grid.
  } else {
    # --- Error handling for higher dimensions ---
    stop("plot.BKP() only supports data where the dimensionality of X is 1 or 2.")
  }
}
