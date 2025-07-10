#' @name plot
#'
#' @title Plot a Beta Kernel Process (BKP) Model
#'
#' @description Visualizes the fitted Beta Kernel Process (BKP) model according
#'   to the dimensionality of the covariate matrix `X`. For 1-dimensional input,
#'   the function plots the estimated probability curve with confidence
#'   intervals and observed data points. For 2-dimensional input, it produces
#'   contour plots for the estimated mean, confidence bounds, and uncertainty.
#'
#' @param x An object of class \code{"BKP"}, typically returned by the
#'   \code{\link{fit.BKP}} function.
#' @param ... Additional arguments passed to generic plot functions (currently
#'   not used, included for S3 method consistency).
#'
#' @details The plotting behavior depends on the dimensionality of the
#'   covariates \code{X} in the fitted BKP model:
#' \itemize{
#'   \item \strong{1D covariates:} Produces a line plot of the predicted probability across the input domain.
#'         A shaded region represents the confidence interval (default 95%), and the observed proportions \code{y/m}
#'         are added as points.
#'   \item \strong{2D covariates:} Generates a 2x2 grid of contour plots showing the predicted mean, upper bound,
#'         lower bound, and the width of the confidence interval (uncertainty).
#' }
#'
#'   If the dimension of \code{X} exceeds 2, the function will stop with an
#'   informative error.
#'
#' @author Jiangyan Zhao, Kunhai Qing, Jin Xu
#'
#' @seealso \code{\link{fit.BKP}} for fitting a BKP model.
#'   \code{\link{predict.BKP}} for generating predictions from a BKP model.
#'
#' @keywords BKP
#'
#' @examples
#' ### 1D
#' set.seed(123)
#' n <- 100
#' Xbounds <- matrix(c(-2,2), nrow=1)
#' x <- seq(-2, 2, length = n)
#' true_pi <- (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#' df <- data.frame(x = x, y = y, m = m)
#' xx = matrix(seq(-2, 2, length = 100), ncol=1) #new data points
#' model <- fit.BKP(df, Xbounds=Xbounds)
#' head(predict(model,xx))
#' plot(model)
#'
#' ### 2D
#' set.seed(123)
#' n <- 100
#' f <- function(X) {
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
#'   return(f) }
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' library(tgp)
#' x <- lhs(n = n, rect = Xbounds)
#' true_pi <- pnorm(f(x))
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#' df <- data.frame(x = x, y = y, m = m)
#' xx1 <- seq(Xbounds[1,1], Xbounds[1,2], length.out = 100)
#' xx2 <- seq(Xbounds[2,1], Xbounds[2,2], length.out = 100)
#' xx <- expand.grid(xx1 = xx1, xx2 = xx2)
#' #plot the true probability
#' true_pi <- matrix(pnorm(f(xx)), nrow = length(xx1), ncol = length(xx2))
#' image(xx1, xx2, true_pi, xlab ="X1", ylab ="X2",
#'                 main = "True Probability",
#'                 col = hcl.colors(100, "viridis"))
#' contour(xx1, xx2, true_pi, add = TRUE, col = "black")
#' model <- fit.BKP(df)
#' head(predict(model,xx))
#' plot(model)
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
    pred <- predict.BKP(BKPmodel, Xnew, CI_size = 0.05)

    # Initialize the plot with the estimated probability curve.
    plot(pred$x1, pred$mean,
         type = "l", col = "blue", lwd = 2,
         xlab = "x (Input Variable)", ylab = "Probability",
         main = "Estimated Probability",
         xlim = Xbounds,
         ylim = c(max(0, min(pred$mean)-0.1),
                  min(1, max(pred$mean)+0.1)))

    # Add a shaded confidence interval band using polygon.
    polygon(c(pred$x1, rev(pred$x1)),
            c(pred$lower, rev(pred$upper)),
            col = "lightgrey", border = NA)
    lines(pred$x1, pred$mean, col = "blue", lwd = 2)

    # Overlay observed proportions (y/m) as points.
    points(X, y / m, pch = 20, col = "red")

    # Add a legend to explain plot elements.
    legend("topright",
           legend = c("Estimated Probability", "95% Confidence Interval", "Observed Proportions"),
           col = c("blue", "lightgrey", "red"), bty = "n",
           lwd = c(2, 8, NA), pch = c(NA, NA, 20), lty = c(1, 1, NA))
  } else if (d == 2){
    #----- Plotting for 2-dimensional covariate data (d == 2) -----#
    # Generate 2D prediction grid
    x1 <- seq(Xbounds[1, 1], Xbounds[1, 2], length.out = 100)
    x2 <- seq(Xbounds[2, 1], Xbounds[2, 2], length.out = 100)
    grid <- expand.grid(x1 = x1, x2 = x2)
    pred <- predict.BKP(BKPmodel, as.matrix(grid), CI_size = 0.05)

    # Convert to data frame for levelplot
    df <- data.frame(x1 = pred$x1, x2 = pred$x2,
                     Mean = pred$mean,
                     Upper = pred$upper,
                     Lower = pred$lower,
                     Variance = pred$variance)
                     # Width = pred$upper - pred$lower)

    # Helper to build each plot
    plot_fun <- function(var, title, pal = "plasma", ...) {
      levelplot(
        as.formula(paste(var, "~ x1 * x2")),
        data = df,
        col.regions = hcl.colors(100, palette = pal),
        main = title,
        xlab = "X1", ylab = "X2",
        contour = TRUE,
        colorkey = TRUE,
        cuts = 15,
        pretty = TRUE,
        scales = list(draw = TRUE, tck = c(1, 0)),
        panel = function(...) {
          panel.levelplot(...)
          panel.contourplot(..., col = "black", lwd = 0.5)
        }
      )
    }

    # Create 4 plots
    p1 <- plot_fun("Mean", "Predictive Mean")
    p2 <- plot_fun("Upper", "95% CI Upper")
    p3 <- plot_fun("Variance", "Predictive Variance")
    # p3 <- plot_fun("Width", "CI Width")
    p4 <- plot_fun("Lower", "95% CI Lower")

    # Arrange into 2×2 layout
    grid.arrange(p1, p2, p3, p4, ncol = 2)
  } else {
    # --- Error handling for higher dimensions ---
    stop("plot.BKP() only supports data where the dimensionality of X is 1 or 2.")
  }
}
