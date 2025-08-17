#' @name plot
#'
#' @title Plot Fitted BKP or DKP Models
#'
#' @description Visualizes fitted \code{BKP} or \code{DKP} models depending on
#'   the input dimensionality. For 1-dimensional inputs, it displays predicted
#'   class probabilities with credible intervals and observed data. For
#'   2-dimensional inputs, it generates contour plots of posterior summaries.
#'
#' @param x An object of class \code{"BKP"} or \code{"DKP"}, typically returned
#'   by \code{\link{fit.BKP}} or \code{\link{fit.DKP}}.
#' @param only_mean Logical. If \code{TRUE}, only the predicted mean surface is
#'   plotted for 2D inputs (only applies to \code{BKP} models). Default is
#'   \code{FALSE}.
#' @param n_grid Integer. Number of grid points used in each dimension to
#'   construct the prediction grid. A larger value produces a smoother and more
#'   detailed decision boundary, but increases computational cost.
#'   Default is \code{80}.
#' @param ... Additional arguments passed to internal plotting routines
#'   (currently unused).
#'
#' @return This function does not return a value. It is called for its side
#'   effects, producing plots that visualize the model predictions and
#'   uncertainty.
#'
#' @details The plotting behavior depends on the dimensionality of the input
#'   covariates:
#'
#' \itemize{
#'   \item \strong{1D inputs:}
#'     \itemize{
#'       \item For \code{BKP} (binary/binomial data), plots the posterior mean curve with a 95% credible band, overlaid with the observed proportions (\eqn{y/m}).
#'       \item For \code{DKP} (categorical/multinomial data), plots one curve per class, each with a shaded credible interval and the observed class frequencies.
#'       \item For classification tasks, an optional curve of the maximum posterior class probability can be displayed to visualize decision confidence.
#'     }
#'
#'   \item \strong{2D inputs:}
#'     \itemize{
#'       \item For both models, produces either:
#'         \enumerate{
#'           \item A predictive mean surface (optionally maximum posterior probability for classification), or
#'           \item A 2-by-2 panel of contour plots showing: predictive mean, predictive 97.5th percentile (upper bound of 95% credible interval), predictive variance, and predictive 2.5th percentile (lower bound).
#'         }
#'       \item For DKP, these surfaces are generated separately for each class.
#'     }
#' }
#'
#'   For input dimensions greater than two, the function terminates with an
#'   error message.
#'
#' @seealso \code{\link{fit.BKP}}, \code{\link{predict.BKP}},
#'   \code{\link{fit.DKP}}, \code{\link{predict.DKP}}
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
#' model1 <- fit.BKP(X, y, m, Xbounds=Xbounds)
#'
#' # Plot results
#' plot(model1)
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
#'
#' # Plot results
#' plot(model2, n_grid = 50)
#'
#' @export
#' @method plot BKP

plot.BKP <- function(x, only_mean = FALSE, n_grid = 80, ...){
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
    Xnew <- matrix(seq(Xbounds[1], Xbounds[2], length.out = 10 * n_grid), ncol = 1)

    # Get the prediction for the new X values.
    prediction <- predict.BKP(BKPmodel, Xnew, ...)
    is_classification <- !is.null(prediction$class)

    # Initialize the plot with the estimated probability curve.
    plot(Xnew, prediction$mean,
         type = "l", col = "blue", lwd = 2,
         xlab = "x", ylab = "Probability",
         main = "Estimated Probability",
         xlim = Xbounds,
         ylim = c(min(prediction$lower) * 0.9,
                  min(1, max(prediction$upper) * 1.1)))

    # Add a shaded credible interval band using polygon.
    polygon(c(Xnew, rev(Xnew)),
            c(prediction$lower, rev(prediction$upper)),
            col = "lightgrey", border = NA)
    lines(Xnew, prediction$mean, col = "blue", lwd = 2)

    # Overlay observed proportions (y/m) as points.
    points(X, y / m, pch = 20, col = "red")

    if(is_classification){
      abline(h = prediction$threshold, lty = 2, lwd = 1.2)
      text(x = Xbounds[2],
           y = prediction$threshold + 0.02,
           labels = "threshold",
           adj = c(1, 0.5),
           cex = 0.9,
           col = "black")
    }else{
      # Add a legend to explain plot elements.
      legend("topleft",
             legend = c("Estimated Probability",
                        paste0(prediction$CI_level * 100, "% Credible Interval"),
                        "Observed Proportions"),
             col = c("blue", "lightgrey", "red"), bty = "n",
             lwd = c(2, 8, NA), pch = c(NA, NA, 20), lty = c(1, 1, NA))
    }

  } else if (d == 2){
    #----- Plotting for 2-dimensional covariate data (d == 2) -----#
    # Generate 2D prediction grid
    x1 <- seq(Xbounds[1, 1], Xbounds[1, 2], length.out = n_grid)
    x2 <- seq(Xbounds[2, 1], Xbounds[2, 2], length.out = n_grid)
    grid <- expand.grid(x1 = x1, x2 = x2)
    prediction <- predict.BKP(BKPmodel, as.matrix(grid), ...)
    is_classification <- !is.null(prediction$class)

    # Convert to data frame for levelplot
    df <- data.frame(x1 = grid$x1, x2 = grid$x2,
                     Mean = prediction$mean,
                     Upper = prediction$upper,
                     Lower = prediction$lower,
                     Variance = prediction$variance)
                     # Width = prediction$upper - prediction$lower)

    if (only_mean) {
      # Only plot the predicted mean graphs
      if(is_classification){
        p1 <- my_2D_plot_fun("Mean", "Predicted Class Probability (Predictive Mean)", df, X = X, y = y)
      }else{
        p1 <- my_2D_plot_fun("Mean", "Predictive Mean", df)
      }
      print(p1)
    } else {
      # Create 2 or 4 plots
      if(is_classification){
        p1 <- my_2D_plot_fun("Mean", "Predictive Mean", df, X = X, y = y)
        p3 <- my_2D_plot_fun("Variance", "Predictive Variance", df, X = X, y = y)
      }else{
        p1 <- my_2D_plot_fun("Mean", "Predictive Mean", df)
        p3 <- my_2D_plot_fun("Variance", "Predictive Variance", df)
      }
      if(is_classification){
        # Arrange into 1×2 layout
        grid.arrange(p1, p3, ncol = 2)
      }else{
        p2 <- my_2D_plot_fun("Upper", paste0(prediction$CI_level * 100, "% CI Upper"), df)
        p4 <- my_2D_plot_fun("Lower", paste0(prediction$CI_level * 100, "% CI Lower"), df)
        # Arrange into 2×2 layout
        grid.arrange(p1, p2, p3, p4, ncol = 2)
      }
    }
  } else {
    # --- Error handling for higher dimensions ---
    stop("plot.BKP() only supports data where the dimensionality of X is 1 or 2.")
  }
}
