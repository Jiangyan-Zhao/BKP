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

plot.DKP <- function(x, only_mean = FALSE, n_grid = 80, dims = NULL,
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

    if (engine == "ggplot") {
      class_labels <- factor(rep(seq_len(q), each = nrow(Xnew)), levels = seq_len(q))
      mean_vals <- as.vector(prediction$mean)
      lower_vals <- as.vector(prediction$lower)
      upper_vals <- as.vector(prediction$upper)
      x_vals <- rep(as.numeric(Xnew), times = q)

      obs_list <- lapply(seq_len(q), function(j) {
        if (is_classification) {
          as.integer(apply(Y, 1, which.max) == j)
        } else {
          Y[, j] / rowSums(Y)
        }
      })
      obs_vals <- unlist(obs_list)
      x_obs <- rep(as.numeric(X_sub), times = q)

      pred_df <- data.frame(
        x = x_vals,
        mean = mean_vals,
        lower = lower_vals,
        upper = upper_vals,
        class = class_labels
      )

      obs_df <- data.frame(
        x = x_obs,
        obs = obs_vals,
        class = factor(rep(seq_len(q), each = nrow(X_sub)), levels = seq_len(q))
      )

      plot_df <- data.frame(
        x = c(x_vals, x_obs),
        mean = c(mean_vals, rep(NA_real_, length(obs_vals))),
        lower = c(lower_vals, rep(NA_real_, length(obs_vals))),
        upper = c(upper_vals, rep(NA_real_, length(obs_vals))),
        class = factor(c(as.character(class_labels), rep(seq_len(q), each = nrow(X_sub))),
                       levels = seq_len(q)),
        obs = c(rep(NA_real_, length(mean_vals)), obs_vals)
      )

      y_max <- if (is_classification) 1 else min(1, max(plot_df$upper, na.rm = TRUE) * 1.1)
      y_min <- if (is_classification) 0 else min(plot_df$lower, na.rm = TRUE) * 0.9

      ci_label <- paste0(prediction$CI_level * 100, "% CI")

      p <- ggplot(plot_df, aes(x = .data$x)) +
        geom_ribbon(data = pred_df,
                    aes(ymin = .data$lower, ymax = .data$upper, fill = ci_label),
                    alpha = 0.3,
                    show.legend = TRUE) +
        geom_line(data = pred_df,
                  aes(y = .data$mean, color = "Estimated Probability"),
                  linewidth = 1,
                  show.legend = TRUE) +
        geom_point(data = obs_df,
                   aes(y = .data$obs, color = "Observed"),
                   size = 1.5,
                   alpha = 0.85,
                   show.legend = TRUE) +
        facet_wrap(~class, ncol = 2, scales = "fixed") +
        coord_cartesian(ylim = c(y_min, y_max)) +
        scale_color_manual(
          values = c("Estimated Probability" = "blue", "Observed" = "red"),
          name = NULL
        ) +
        scale_fill_manual(values = setNames("grey70", ci_label), name = NULL) +
        guides(color = guide_legend(order = 1), fill = guide_legend(order = 2)) +
        labs(
          title = "Estimated Probability by Class",
          x = if (is.null(dims)) "x" else paste0("x", dims[1]),
          y = "Probability"
        ) +
        theme_minimal()
      print(p)
    } else {
      old_par <- par(mfrow = c(2, 2))
      # on.exit(par(old_par))  # Restore par on exit

      # --- First panel: all mean curves together ---
      if(is_classification){
        cols <- rainbow(q)
        plot(NA,
             xlim = Xbounds[dims, ],
             ylim = c(-0.1, 1.1),
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
             xlim = Xbounds[dims, ],
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

      if (engine == "ggplot") {
        p1 <- my_2D_plot_fun_class_ggplot("class", "Predicted Classes", df, X_sub, Y, dims = dims)
        p2 <- my_2D_plot_fun_class_ggplot("max_prob", "Maximum Predicted Probability", df, X_sub, Y, classification = FALSE, dims = dims)
      } else {
        p1 <- my_2D_plot_fun_class("class", "Predicted Classes", df, X_sub, Y, dims= dims)
        p2 <- my_2D_plot_fun_class("max_prob", "Maximum Predicted Probability", df, X_sub, Y, classification = FALSE, dims= dims)
      }
      gridExtra::grid.arrange(p1, p2, ncol = 2)
    }else{
      for (j in 1:q) {
        df <- data.frame(x1 = grid$x1, x2 = grid$x2,
                         Mean = prediction$mean[, j],
                         Upper = prediction$upper[, j],
                         Lower = prediction$lower[, j],
                         Variance = prediction$variance[, j])

        if (only_mean) {
          # Only plot the predicted mean graphs
          p1 <- if (engine == "ggplot") {
            my_2D_plot_fun_ggplot("Mean", "Predictive Mean", df, dims = dims)
          } else {
            my_2D_plot_fun("Mean", "Predictive Mean", df, dims= dims)
          }
          print(p1)
        } else {
          # Create 4 plots
          if (engine == "ggplot") {
            p1 <- my_2D_plot_fun_ggplot("Mean", "Predictive Mean", df, dims = dims)
            p2 <- my_2D_plot_fun_ggplot("Upper", paste0(prediction$CI_level * 100, "% CI Upper"), df, dims = dims)
            p3 <- my_2D_plot_fun_ggplot("Variance", "Predictive Variance", df, dims = dims)
            p4 <- my_2D_plot_fun_ggplot("Lower", paste0(prediction$CI_level * 100, "% CI Lower"), df, dims = dims)
          } else {
            p1 <- my_2D_plot_fun("Mean", "Predictive Mean", df, dims= dims)
            p2 <- my_2D_plot_fun("Upper", paste0(prediction$CI_level * 100, "% CI Upper"), df, dims= dims)
            p3 <- my_2D_plot_fun("Variance", "Predictive Variance", df, dims= dims)
            p4 <- my_2D_plot_fun("Lower", paste0(prediction$CI_level * 100, "% CI Lower"), df, dims= dims)
          }
          # Arrange into 2×2 layout
          gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2,
                       top = textGrob(paste0("Estimated Probability (Class ", j, ")"),
                                      gp = gpar(fontface = "bold", fontsize = 16)))
        }
      }
    }
  }

  # par(old_par)
}
