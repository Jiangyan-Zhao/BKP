#' @name plot
#'
#' @title Plot a Dirichlet Kernel Process (DKP) Model
#'
#' @description Visualizes the fitted Dirichlet Kernel Process (DKP) model according
#'   to the dimensionality of the covariate matrix `X`. For 1-dimensional input,
#'   the function plots the estimated probability curves with confidence
#'   intervals and observed data points. For 2-dimensional input, it produces
#'   contour plots for the estimated mean, confidence bounds, and uncertainty.
#'
#' @param x An object of class \code{"DKP"}, typically returned by the
#'   \code{\link{fit.DKP}} function.
#' @param ... Additional arguments passed to generic plot functions (currently
#'   not used, included for S3 method consistency).
#'
#' @details The plotting behavior depends on the dimensionality of the
#'   covariates \code{X} in the fitted DKP model:
#' \itemize{
#'   \item \strong{1D covariates:} Produces line plots of the predicted probability across the input domain.
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
#' @seealso \code{\link{fit.DKP}} for fitting a DKP model.
#'   \code{\link{predict.DKP}} for generating predictions from a DKP model.
#'
#' @keywords DKP
#'
#' @examples
#' ### 1D
#' set.seed(123)
#' n <- 30
#' Xbounds <- matrix(c(-2,2), nrow=1)
#' x <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' true_pi <- matrix(c(true_pi/2,true_pi/2,1-true_pi),nrow = n, byrow = FALSE)
#' m <- sample(100, n, replace = TRUE)
#' Y <- matrix(0, nrow = n, ncol = 3)
#' for (i in 1:n) {
#'   Y[i, ] <- rmultinom(n=1, size=m[i], prob=true_pi[i, ])
#' }
#' DKPmodel <- fit.DKP(x, Y, Xbounds = Xbounds)
#' plot(DKPmodel)
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
#' x <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- pnorm(f(x))
#' true_pi <- matrix(c(true_pi/2,true_pi/2,1-true_pi),nrow = n, byrow = FALSE)
#' m <- sample(100, n, replace = TRUE)
#' Y <- matrix(0, nrow = n, ncol = 3)
#' for (i in 1:n) {
#'   Y[i, ] <- rmultinom(n=1, size=m[i], prob=true_pi[i, ])
#' }
#' DKPmodel <- fit.DKP(x, Y, Xbounds = Xbounds)
#' plot(DKPmodel)
#'
#' @export
#' @method plot DKP

plot.DKP <- function(x, ...){
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
    Xnew <- matrix(seq(Xbounds[1], Xbounds[2], length.out = 1000), ncol = 1)

    # Get the prediction for the new X values.
    prediction <- predict.DKP(DKPmodel, Xnew)

    for (j in 1:q) {
      mean_j <- prediction$mean[, j]
      lower_j <- prediction$lower[, j]
      upper_j <- prediction$upper[, j]
      # Initialize the plot with the estimated probability curve.
      plot(Xnew, mean_j,
           type = "l", col = "blue", lwd = 2,
           xlab = "x (Input Variable)", ylab = "Probability",
           main = paste0("Estimated Probability for Y Label ",j),
           xlim = Xbounds,
           ylim = c(max(0, min(mean_j)-0.1),
                    min(1, max(mean_j)+0.1)))

      # Add a shaded confidence interval band using polygon.
      polygon(c(Xnew, rev(Xnew)),
              c(lower_j, rev(upper_j)),
              col = "lightgrey", border = NA)
      lines(Xnew, mean_j, col = "blue", lwd = 2)

      # Overlay observed proportions (y/m) as points.
      points(X, Y[,j] / rowSums(Y), pch = 20, col = "red")

      # Add a legend to explain plot elements.
      legend("topright",
             legend = c("Estimated Probability", "95% Confidence Interval", "Observed Proportions"),
             col = c("blue", "lightgrey", "red"), bty = "n",
             lwd = c(2, 8, NA), pch = c(NA, NA, 20), lty = c(1, 1, NA))
    }



  } else if (d == 2){
    #----- Plotting for 2-dimensional covariate data (d == 2) -----#
    # Generate 2D prediction grid
    x1 <- seq(Xbounds[1, 1], Xbounds[1, 2], length.out = 100)
    x2 <- seq(Xbounds[2, 1], Xbounds[2, 2], length.out = 100)
    grid <- expand.grid(x1 = x1, x2 = x2)
    prediction <- predict.DKP(DKPmodel, as.matrix(grid))

    for (j in 1:q) {
      df <- data.frame(x1 = grid$x1, x2 = grid$x2,
                       Mean = prediction$mean[, j],
                       Upper = prediction$upper[, j],
                       Lower = prediction$lower[, j],
                       Variance = prediction$variance[, j])

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
      p4 <- plot_fun("Lower", "95% CI Lower")

      # Arrange into 2×2 layout
      grid.arrange(p1, p2, p3, p4, ncol = 2,
                   top = textGrob(paste0("Estimated Probability for Y Label ",j),
                                  gp = gpar(fontface = "bold", fontsize = 16)))
    }
  } else {
    # --- Error handling for higher dimensions ---
    stop("plot.DKP() only supports data where the dimensionality of X is 1 or 2.")
  }
}
