#' @name plot
#'
#' @title Plot Fitted BKP or DKP Models
#'
#' @description Visualizes fitted \code{BKP} (binary) or \code{DKP}
#'   (multi-class) models according to the input dimensionality. For 1D inputs,
#'   it shows predicted class probabilities with credible intervals and observed
#'   data. For 2D inputs, it generates contour plots of posterior summaries. For
#'   higher-dimensional inputs, users must specify which dimensions to plot.
#'
#' @param x An object of class \code{"BKP"} or \code{"DKP"}, typically returned
#'   by \code{\link{fit_BKP}} or \code{\link{fit_DKP}}.
#' @param only_mean Logical. If \code{TRUE}, only the predicted mean surface is
#'   plotted for 2D inputs (applies to both BKP and DKP models for mean
#'   visualization). Default is \code{FALSE}.
#' @param n_grid Positive integer specifying the number of grid points per
#'   dimension for constructing the prediction grid. Larger values produce
#'   smoother and more detailed surfaces, but increase computation time. Default
#'   is \code{80}.
#' @param dims Integer vector indicating which input dimensions to plot. Must
#'   have length 1 (for 1D) or 2 (for 2D). If \code{NULL} (default), all
#'   dimensions are used when their number is <= 2.
#' @param engine Character string specifying plotting backend for 2D plots.
#'   Either \code{"base"} (default, lattice-based output) or
#'   \code{"ggplot"}.
#' @param ... Additional arguments passed to internal plotting routines
#'   (currently unused).
#'
#' @return This function is called for its side effects and does not return a
#'   value. It produces plots visualizing the predicted probabilities, credible
#'   intervals, and posterior summaries.
#'
#' @details
#' The plotting behavior depends on the dimensionality of the input covariates:
#'
#' \itemize{
#'   \item \strong{1D inputs:}
#'     \itemize{
#'       \item For \code{BKP} (binary/binomial data), the function plots the posterior mean curve with a shaded 95% credible interval, overlaid with the observed proportions (\eqn{y/m}).
#'       \item For \code{DKP} (categorical/multinomial data), it plots one curve per class, each with a shaded credible interval and the observed class frequencies.
#'       \item For classification tasks, an optional curve of the maximum posterior class probability can be displayed to visualize decision confidence.
#'     }
#'
#'   \item \strong{2D inputs:}
#'     \itemize{
#'       \item For both BKP and DKP models, the function generates contour plots over a 2D prediction grid.
#'       \item Users can choose to plot only the predictive mean surface (\code{only_mean = TRUE}) or a set of four summary plots (\code{only_mean = FALSE}):
#'         \enumerate{
#'           \item Predictive mean
#'           \item 97.5th percentile (upper bound of 95% credible interval)
#'           \item Predictive variance
#'           \item 2.5th percentile (lower bound of 95% credible interval)
#'         }
#'       \item For DKP, these surfaces are generated separately for each class.
#'       \item For classification tasks, predictive class probabilities can also be visualized as the maximum posterior probability surface.
#'     }
#'
#'   \item \strong{Input dimensions greater than 2:}
#'     \itemize{
#'       \item The function does not automatically support visualization and will terminate with an error.
#'       \item Users must specify which dimensions to visualize via the \code{dims} argument (length 1 or 2).
#'     }
#' }
#'
#' @seealso \code{\link{fit_BKP}} and \code{\link{fit_DKP}} for fitting BKP and
#'   DKP models, respectively; \code{\link{predict.BKP}} and
#'   \code{\link{predict.DKP}} for generating predictions from fitted BKP and
#'   DKP models.

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
#' model2 <- fit_BKP(X, y, m, Xbounds=Xbounds)
#'
#' # Plot results
#' plot(model2, n_grid = 50)
#'
#' @export
#' @method plot BKP

plot.BKP <- function(x, only_mean = FALSE, n_grid = 80, dims = NULL,
                     engine = c("base", "ggplot"), ...){
  # ---------------- Argument Checking ----------------
  if (!is.logical(only_mean) || length(only_mean) != 1) {
    stop("`only_mean` must be a single logical value (TRUE or FALSE).")
  }

  if (!is.numeric(n_grid) || length(n_grid) != 1 || n_grid <= 0) {
    stop("'n_grid' must be a positive integer.")
  }
  n_grid <- as.integer(n_grid)

  engine <- match.arg(engine)

  # Extract necessary components from the BKP model object.
  X <- x$X # Covariate matrix.
  y <- x$y # Number of successes.
  m <- x$m # Number of trials.
  Xbounds <- x$Xbounds

  d <- ncol(X)    # Dimensionality.

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

  if (length(dims) == 1){
    #----- Plotting for 1-dimensional covariate data (d == 1) -----#
    # Generate new X values for a smooth prediction curve.
    Xnew <- matrix(seq(Xbounds[dims, 1], Xbounds[dims, 2], length.out = 10 * n_grid), ncol = 1)

    # Get the prediction for the new X values.
    Xnew_full <- lhs(nrow(Xnew), Xbounds)
    Xnew_full[, dims] <- Xnew
    prediction <- predict.BKP(x, Xnew_full, ...)

    # Determine whether it is a classification problem
    is_classification <- !is.null(prediction$class)

    if (engine == "ggplot") {
      plot_df <- data.frame(
        x = as.numeric(Xnew),
        mean = prediction$mean,
        lower = prediction$lower,
        upper = prediction$upper
      )
      obs_df <- data.frame(
        x = as.numeric(X_sub),
        obs_y = y / m
      )
      lbl_line <- "Estimated Probability"
      lbl_ci   <- paste0(prediction$CI_level * 100, "% Credible Interval")
      lbl_pts  <- "Observed Proportions"

      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x)) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
          fill = "grey70", alpha = 0.4
        ) +
        ggplot2::geom_line(
          ggplot2::aes(y = .data$mean, color = lbl_ci),
          alpha = 0
        ) +
        ggplot2::geom_line(
          ggplot2::aes(y = .data$mean, color = lbl_line),
          linewidth = 1
        ) +
        ggplot2::geom_point(
          data = obs_df,
          ggplot2::aes(x = .data$x, y = .data$obs_y, color = lbl_pts),
          size = 2
        ) +
        ggplot2::scale_color_manual(
          name = NULL,
          values = stats::setNames(c("blue", "grey70", "red"), c(lbl_line, lbl_ci, lbl_pts)),
          breaks = c(lbl_line, lbl_ci, lbl_pts)
        ) +
        ggplot2::guides(
          color = ggplot2::guide_legend(
            override.aes = list(
              shape     = c(NA, NA, 16),
              linetype  = c(1, 1, 0),
              linewidth = c(1, 5, 0),
              alpha     = c(1, 0.5, 1)
            )
          )
        ) +
        ggplot2::labs(
          title = "Estimated Probability",
          x = ifelse(d > 1, paste0("x", dims), "x"),
          y = "Probability"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(colour = "black", fill=NA, linewidth=1),
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 15),
          axis.title = ggplot2::element_text(size = 13),
          axis.text  = ggplot2::element_text(size = 11, color = "black"),
          legend.position = c(0.02, 0.98),
          legend.justification = c(0, 1),
          legend.background = ggplot2::element_blank(),
          legend.key = ggplot2::element_blank(),
          legend.text = ggplot2::element_text(size = 12),
          legend.key.width = ggplot2::unit(2, "line")
        )

      if (is_classification) {
        p <- p +
          geom_hline(yintercept = prediction$threshold, linetype = "dashed") +
          annotate("text",
                   x = max(plot_df$x), y = prediction$threshold + 0.02,
                   label = "threshold", hjust = 1, vjust = 0.5)
      }
      print(p)
    } else {
      # Initialize the plot with the estimated probability curve.
      plot(Xnew, prediction$mean,
           type = "l", col = "blue", lwd = 2,
           xlab = ifelse(d > 1, paste0("x", dims), "x"),
           ylab = "Probability",
           main = "Estimated Probability",
           xlim = Xbounds[dims, ],
           ylim = c(min(prediction$lower) * 0.9,
                    min(1, max(prediction$upper) * 1.1)))

      # Add a shaded credible interval band using polygon.
      polygon(c(Xnew, rev(Xnew)),
              c(prediction$lower, rev(prediction$upper)),
              col = "lightgrey", border = NA)
      lines(Xnew, prediction$mean, col = "blue", lwd = 2)

      # Overlay observed proportions (y/m) as points.
      points(X_sub, y / m, pch = 20, col = "red")

      if(is_classification){
        abline(h = prediction$threshold, lty = 2, lwd = 1.2)
        text(x = Xbounds[dims, 2],
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
    prediction <- predict.BKP(x, Xnew_full, ...)

    # Determine whether it is a classification problem
    is_classification <- !is.null(prediction$class)

    # Convert to data frame for levelplot
    df <- data.frame(x1 = grid$x1, x2 = grid$x2,
                     Mean = prediction$mean,
                     Upper = prediction$upper,
                     Lower = prediction$lower,
                     Variance = prediction$variance)

    if (only_mean) {
      # Only plot the predicted mean graphs
      if(is_classification){
        p1 <- if (engine == "ggplot") {
          my_2D_plot_fun_ggplot("Mean", "Predicted Class Probability (Predictive Mean)", df, X = X_sub, y = y, dims = dims)
        } else {
          my_2D_plot_fun("Mean", "Predicted Class Probability (Predictive Mean)", df, X = X_sub, y = y, dims = dims)
        }
      }else{
        p1 <- if (engine == "ggplot") {
          my_2D_plot_fun_ggplot("Mean", "Predictive Mean", df, dims = dims)
        } else {
          my_2D_plot_fun("Mean", "Predictive Mean", df, dims = dims)
        }
      }
      print(p1)
    } else {
      # Create 2 or 4 plots
      if(is_classification){
        if (engine == "ggplot") {
          p1 <- my_2D_plot_fun_ggplot("Mean", "Predictive Mean", df, X = X_sub, y = y, dims = dims)
          p3 <- my_2D_plot_fun_ggplot("Variance", "Predictive Variance", df, X = X_sub, y = y, dims = dims)
        } else {
          p1 <- my_2D_plot_fun("Mean", "Predictive Mean", df, X = X_sub, y = y, dims= dims)
          p3 <- my_2D_plot_fun("Variance", "Predictive Variance", df, X = X_sub, y = y, dims= dims)
        }
      }else{
        if (engine == "ggplot") {
          p1 <- my_2D_plot_fun_ggplot("Mean", "Predictive Mean", df, dims = dims)
          p3 <- my_2D_plot_fun_ggplot("Variance", "Predictive Variance", df, dims = dims)
        } else {
          p1 <- my_2D_plot_fun("Mean", "Predictive Mean", df, dims= dims)
          p3 <- my_2D_plot_fun("Variance", "Predictive Variance", df, dims= dims)
        }
      }
      if(is_classification){
        # Arrange into 1×2 layout
        gridExtra::grid.arrange(p1, p3, ncol = 2)
      }else{
        if (engine == "ggplot") {
          p2 <- my_2D_plot_fun_ggplot("Upper", paste0(prediction$CI_level * 100, "% CI Upper"), df, dims = dims)
          p4 <- my_2D_plot_fun_ggplot("Lower", paste0(prediction$CI_level * 100, "% CI Lower"), df, dims = dims)
        } else {
          p2 <- my_2D_plot_fun("Upper", paste0(prediction$CI_level * 100, "% CI Upper"), df, dims= dims)
          p4 <- my_2D_plot_fun("Lower", paste0(prediction$CI_level * 100, "% CI Lower"), df, dims= dims)
        }
        # Arrange into 2×2 layout
        gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
      }
    }
  }
}
