#' @title Plot a Beta Kernel Process (BKP) model
#'
#' @description
#' This function visualizes the fitted BKP model based on the dimensionality of the covariate data (X).
#' For 1-dimensional covariate data, it plots the estimated probability curve along with a
#' confidence interval and overlays the observed data points.
#' For 2-dimensional covariate data, it generates a set of contour plots displaying the
#' estimated probability, upper confidence bound, lower confidence bound, and confidence interval width.
#'
#' @param x An object of class "bkp", typically returned by the `fit.bkp()` function.
#' @param ... Additional arguments passed to generic plot functions (currently not used
#'   by this specific method but included for S3 generic consistency).
#'
#' @return This function is called for its side effect of generating plots. It returns `NULL` invisibly.
#'
#' @details
#' The plot behavior adapts to the dimension of the covariates (`X`) in the `bkp` model:
#' \itemize{
#'   \item **1-Dimensional X:** A line plot shows the predicted probability (`piHat`) against X.
#'     A shaded area represents the confidence interval (default 95%), and observed proportions (`y/m`)
#'     are added as points.
#'   \item **2-Dimensional X:** A 2x2 grid of contour plots is generated. Each plot visualizes
#'     a different aspect across the 2D covariate space: estimated probability, upper confidence
#'     bound, lower confidence bound, and the width of the confidence interval. Color palettes
#'     are used to represent values.
#' }
#' If `X` has more than 2 dimensions, the function will stop with an error.
#'
#'
#' @seealso
#' \code{\link{fit.bkp}} to fit a BKP model.
#' \code{\link{predict.bkp}} to make predictions from a BKP model.
#'
#' @export
#' @method plot bkp
#' @importFrom graphics plot polygon points image layout legend contour
#' @importFrom grDevices rgb hcl.colors
#'
#' @examples
#' ### 1D
#' set.seed(123)
#' n <- 100
#' x <- seq(-2, 2, length = n)
#' true_pi <- (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#' df <- data.frame(x = x, y = y, m = m)
#' xx <- seq(-2, 2, length = 100) #new data points
#' model <- fit.bkp(df)
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
#' x <- lhs::randomLHS(n = n, k = 2)
#' true_pi <- pnorm(f(x))
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#' df <- data.frame(x = x, y = y, m = m)
#' xx1 <- seq(0, 1, length.out = 100)
#' xx2 <- seq(0, 1, length.out = 100)
#' xx <- expand.grid(xx1 = xx1, xx2 = xx2)
#' model <- fit.bkp(df)
#' plot(model)
#'
#' #plot the true probability
#' true_pi <- pnorm(f(xx))
#' pred_prob_matrix <- matrix(true_pi,
#'                            nrow = length(xx1),
#'                            ncol = length(xx2))
#' graphics::image(xx1, xx2, pred_prob_matrix,
#'                 xlab = "X1", # Use X column name if available.
#'                 ylab = "X2", # Use X column name if available.
#'                 main = "True Probability",
#'                 col = grDevices::hcl.colors(100, "viridis"))
plot.bkp <- function(
    x,
    ...
){
  bkp.model <- x # Assign the input 'x' (the bkp object) to a more descriptive variable.

  # Check if the provided object is indeed of class "bkp".
  # This ensures the function operates on a valid BKP model output.
  if (!inherits(bkp.model, "bkp")) {
    stop("The object in question is not of class \"bkp\". Please provide a model fitted by 'fit.bkp()'.")
  }

  # Extract necessary components from the BKP model object.
  y <- bkp.model$y # Number of successes.
  m <- bkp.model$m # Number of trials.
  X <- bkp.model$originalX # Original (unnormalized) covariate matrix.
  d <- ncol(X) # Determine the dimensionality of the covariates.

  # --- Plotting for 1-dimensional covariate data (d == 1) ---
  if (d == 1){
    # Generate new X values for a smooth prediction curve.
    Xnew <- seq(min(X), max(X), length.out = 1000)
    # Predict probabilities and confidence intervals for the new X values.
    # Assumes 'predict.bkp' is an existing method.
    prediction <- predict.bkp(bkp.model, Xnew, CI.size = 0.05)

    # Initialize the plot with the estimated probability curve.
    graphics::plot(prediction$X, prediction$piHat,
                   type = "l", col = "blue", lwd = 2, # Line type, color, and width.
                   ylab = "Probability", # Y-axis label.
                   xlab = ifelse(!is.null(colnames(X)[1]), colnames(X)[1], "X"), # Use X column name if available.
                   main = "BKP Estimated Probability", # Main plot title.
                   ylim = c(0, 1)) # Ensure Y-axis spans the full probability range.

    # Add a shaded confidence interval band using polygon.
    graphics::polygon(c(prediction$X, rev(prediction$X)), # X coordinates for polygon.
                      c(prediction$lower, rev(prediction$upper)), # Y coordinates (lower then reversed upper).
                      col = grDevices::rgb(0.6, 0.8, 1, 0.5), # Light blue color with transparency.
                      border = NA) # No border for the polygon.

    # Overlay observed proportions (y/m) as points.
    graphics::points(X, y / m, pch = 16, cex = 0.8, col = "red")

    # Add a legend to explain plot elements.
    graphics::legend("topright",
                     legend = c("Estimated Probability", "95% Confidence Interval", "Observed Proportions"),
                     col = c("blue", grDevices::rgb(0.6, 0.8, 1, 0.5), "red"),
                     lwd = c(2, NA, NA), # Line width for estimated curve.
                     pch = c(NA, 15, 16), # Point/box shapes for CI and observed.
                     lty = c(1, NA, NA), # Line type for estimated curve.
                     bty = "n") # No box around the legend.

    # --- Plotting for 2-dimensional covariate data (d == 2) ---
  } else if (d == 2){
    # Extract individual X dimensions.
    x1 <- X[, 1]
    x2 <- X[, 2]

    # Create a grid of new X values for prediction over the 2D space.
    xx1 <- seq(min(x1), max(x1), length.out = 100)
    xx2 <- seq(min(x2), max(x2), length.out = 100)
    Xnew <- expand.grid(xx1 = xx1, xx2 = xx2) # Create all combinations of xx1 and xx2.

    # Predict probabilities and confidence intervals for the grid.
    prediction <- predict.bkp(bkp.model, Xnew, CI.size = 0.05)

    # Reshape the 1D prediction results into 2D matrices for image plotting.
    pred_prob_matrix <- matrix(prediction$piHat,
                               nrow = length(xx1),
                               ncol = length(xx2))
    pred_upper_matrix <- matrix(prediction$upper,
                                nrow = length(xx1),
                                ncol = length(xx2))
    pred_lower_matrix <- matrix(prediction$lower,
                                nrow = length(xx1),
                                ncol = length(xx2))
    # Calculate the width of the confidence interval.
    pred_upper_lower_matrix <- pred_upper_matrix - pred_lower_matrix

    # Set up a 2x2 plot layout for multiple contour plots.
    mat <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
    graphics::layout(mat)

    # Plot 1: Estimated Probability Contour Map
    graphics::image(xx1, xx2, pred_prob_matrix,
                    xlab = ifelse(!is.null(colnames(X)[1]), colnames(X)[1], "X1"), # Use X column name if available.
                    ylab = ifelse(!is.null(colnames(X)[2]), colnames(X)[2], "X2"), # Use X column name if available.
                    main = "Estimated Probability",
                    col = grDevices::hcl.colors(100, "viridis")) # Use a perceptually uniform color palette.
    graphics::contour(xx1, xx2, pred_prob_matrix, add = TRUE, col = "black") # Overlay black contour lines.

    # Plot 2: Upper Bound Contour Map
    graphics::image(xx1, xx2, pred_upper_matrix,
                    xlab = ifelse(!is.null(colnames(X)[1]), colnames(X)[1], "X1"),
                    ylab = ifelse(!is.null(colnames(X)[2]), colnames(X)[2], "X2"),
                    main = "Upper Bound",
                    col = grDevices::hcl.colors(100, "viridis"))
    graphics::contour(xx1, xx2, pred_upper_matrix, add = TRUE, col = "black")

    # Plot 3: Lower Bound Contour Map
    graphics::image(xx1, xx2, pred_lower_matrix,
                    xlab = ifelse(!is.null(colnames(X)[1]), colnames(X)[1], "X1"),
                    ylab = ifelse(!is.null(colnames(X)[2]), colnames(X)[2], "X2"),
                    main = "Lower Bound",
                    col = grDevices::hcl.colors(100, "viridis"))
    graphics::contour(xx1, xx2, pred_lower_matrix, add = TRUE, col = "black")

    # Plot 4: Confidence Interval Width Contour Map
    graphics::image(xx1, xx2, pred_upper_lower_matrix,
                    xlab = ifelse(!is.null(colnames(X)[1]), colnames(X)[1], "X1"),
                    ylab = ifelse(!is.null(colnames(X)[2]), colnames(X)[2], "X2"),
                    main = "CI Width",
                    col = grDevices::hcl.colors(100, "plasma")) # Use a different palette for width.
    graphics::contour(xx1, xx2, pred_upper_lower_matrix, add = TRUE, col = "black")

    graphics::layout(1) # Reset the plot layout to default (1 plot per page) after generating the grid.

    # --- Error handling for higher dimensions ---
  } else {
    base::stop("plot.bkp() only supports data where the dimensionality of X is 1 or 2.")
  }
  invisible(NULL) # Return NULL invisibly, as the function's primary purpose is plotting.
}
