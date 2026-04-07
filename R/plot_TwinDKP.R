#' @rdname plot
#' @keywords DKP
#'
#' @details
#' \code{plot.TwinDKP()} follows the same workflow as \code{plot.DKP()}, but
#' all predictive quantities are computed via \code{predict.TwinDKP()} on a
#' constructed grid.
#'
#' For 1D: draw class-wise mean and credible interval curves with observed
#' proportions (or class indicators for single-label classification).
#'
#' For 2D: draw predicted class map + max probability for classification, or
#' class-wise mean/upper/variance/lower panels for multinomial probability data.
#'
#' @param l_nums Optional local neighbor count passed to
#'   \code{predict.TwinDKP()}.
#' @param v_nums Optional validation size passed to
#'   \code{predict.TwinDKP()}.
#'
#' @export
#' @method plot TwinDKP
plot.TwinDKP <- function(x, only_mean = FALSE, n_grid = 80, dims = NULL,
                         engine = c("base", "ggplot"),
                         l_nums = NULL, v_nums = NULL, ...) {
  # ---------------- Argument Checking ----------------
  if (!is.logical(only_mean) || length(only_mean) != 1) {
    stop("`only_mean` must be a single logical value (TRUE or FALSE).")
  }

  if (!is.numeric(n_grid) || length(n_grid) != 1 || n_grid <= 0) {
    stop("'n_grid' must be a positive integer.")
  }
  n_grid <- as.integer(n_grid)

  engine <- match.arg(engine)

  # Extract necessary components from the TwinDKP model object.
  X <- x$X
  Y <- x$Y
  Xbounds <- x$Xbounds

  d <- ncol(X)
  q <- ncol(Y)

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

  if (length(dims) == 1) {
    #----- Plotting for 1-dimensional covariate data (d == 1) -----#

    Xnew <- matrix(seq(Xbounds[dims, 1], Xbounds[dims, 2], length.out = 10 * n_grid), ncol = 1)

    Xnew_full <- tgp::lhs(nrow(Xnew), Xbounds)
    Xnew_full[, dims] <- Xnew
    prediction <- predict.TwinDKP(
      x,
      Xnew = Xnew_full,
      l_nums = l_nums,
      v_nums = v_nums,
      ...
    )

    is_classification <- !is.null(prediction$class)

    if (engine == "ggplot") {
      plot_list <- list()

      lbl_line <- "Estimated Probability"
      lbl_ci   <- paste0(prediction$CI_level * 100, "% CI")
      lbl_pts  <- "Observed"

      if (is_classification) {
        all_pred_df <- data.frame(
          x = rep(as.numeric(Xnew), q),
          prob = as.vector(prediction$mean),
          Class = factor(rep(1:q, each = nrow(Xnew)))
        )
        obs_class <- apply(Y, 1, which.max)
        all_obs_df <- data.frame(
          x = as.numeric(X_sub),
          y = rep(-0.05, nrow(X_sub)),
          Class = factor(obs_class, levels = 1:q)
        )

        p_all <- ggplot() +
          geom_line(data = all_pred_df, aes(x = .data$x, y = .data$prob, color = .data$Class), linewidth = 1) +
          geom_point(data = all_obs_df, aes(x = .data$x, y = .data$y, color = .data$Class), size = 2) +
          scale_color_discrete(name = NULL, labels = paste("Class", 1:q)) +
          labs(
            title = "Estimated Mean Curves (All Classes)",
            x = ifelse(d > 1, paste0("x", dims), "x"),
            y = "Probability"
          ) +
          coord_cartesian(ylim = c(-0.1, 1.1)) +
          theme_bw() +
          theme(
            panel.grid = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
            legend.position = "top",
            legend.direction = "horizontal",
            legend.background = element_blank(),
            legend.key = element_blank()
          )
        plot_list[[1]] <- p_all
      }

      for (j in 1:q) {
        mean_j  <- prediction$mean[, j]
        lower_j <- prediction$lower[, j]
        upper_j <- prediction$upper[, j]

        pred_df_j <- data.frame(x = as.numeric(Xnew), mean = mean_j, lower = lower_j, upper = upper_j)

        if (is_classification) {
          obs_j <- as.integer(apply(Y, 1, which.max) == j)
          ylim_j <- c(0, 1)
        } else {
          obs_j <- Y[, j] / rowSums(Y)
          ylim_j <- c(min(lower_j) * 0.9, min(1, max(upper_j) * 1.1))
        }
        obs_df_j <- data.frame(x = as.numeric(X_sub), obs = obs_j)

        p <- ggplot() +
          geom_ribbon(data = pred_df_j, aes(x = .data$x, ymin = .data$lower, ymax = .data$upper), fill = "grey70", alpha = 0.4) +
          geom_line(data = pred_df_j, aes(x = .data$x, y = .data$mean, color = lbl_ci), alpha = 0) +
          geom_line(data = pred_df_j, aes(x = .data$x, y = .data$mean, color = lbl_line), linewidth = 1) +
          geom_point(data = obs_df_j, aes(x = .data$x, y = .data$obs, color = lbl_pts), size = 2) +
          scale_color_manual(name = NULL, values = stats::setNames(c("blue", "grey70", "red"), c(lbl_line, lbl_ci, lbl_pts)), breaks = c(lbl_line, lbl_ci, lbl_pts)) +
          guides(color = guide_legend(override.aes = list(shape = c(NA, NA, 16), linetype = c(1, 1, 0), linewidth = c(1, 5, 0), alpha = c(1, 0.5, 1)))) +
          labs(
            title = paste0("Estimated Probability (Class ", j, ")"),
            x = ifelse(d > 1, paste0("x", dims), "x"),
            y = "Probability"
          ) +
          coord_cartesian(ylim = ylim_j) +
          theme_bw() +
          theme(
            panel.grid = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
            axis.title = element_text(size = 12),
            axis.text  = element_text(size = 10, color = "black")
          )

        if (j == 1) {
          p <- p + theme(
            legend.position = c(0.02, 0.98),
            legend.justification = c(0, 1),
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.text = element_text(size = 11),
            legend.key.width = unit(2, "line")
          )
        } else {
          p <- p + theme(legend.position = "none")
        }

        if (is_classification) {
          plot_list[[j + 1]] <- p
        } else {
          plot_list[[j]] <- p
        }
      }

      do.call(grid.arrange, c(plot_list, ncol = 2))
    } else {
      old_par <- par(mfrow = c(2, 2))

      if (is_classification) {
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

      for (j in 1:q) {
        mean_j  <- prediction$mean[, j]
        lower_j <- prediction$lower[, j]
        upper_j <- prediction$upper[, j]

        if (is_classification) {
          ylim = c(0, 1)
        } else {
          ylim = c(min(lower_j) * 0.9, min(1, max(upper_j) * 1.1))
        }
        plot(Xnew, mean_j,
             type = "l", col = "blue", lwd = 2,
             xlab = ifelse(d > 1, paste0("x", dims), "x"),
             ylab = "Probability",
             main = paste0("Estimated Probability (Class ", j, ")"),
             xlim = Xbounds[dims, ],
             ylim = ylim)

        polygon(c(Xnew, rev(Xnew)),
                c(lower_j, rev(upper_j)),
                col = "lightgrey", border = NA)
        lines(Xnew, mean_j, col = "blue", lwd = 2)

        if (is_classification) {
          obs_j <- as.integer(apply(Y, 1, which.max) == j)
          points(X_sub, obs_j, pch = 20, col = "red")
        } else {
          points(X_sub, Y[, j] / rowSums(Y), pch = 20, col = "red")

          if (j == 1) {
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

      par(old_par)
    }
  } else {
    #----- Plotting for 2-dimensional covariate data (d == 2) -----#
    x1 <- seq(Xbounds[dims[1], 1], Xbounds[dims[1], 2], length.out = n_grid)
    x2 <- seq(Xbounds[dims[2], 1], Xbounds[dims[2], 2], length.out = n_grid)
    grid <- expand.grid(x1 = x1, x2 = x2)

    Xnew_full <- lhs(nrow(grid), Xbounds)
    Xnew_full[, dims] <- as.matrix(grid)
    prediction <- predict.TwinDKP(
      x,
      Xnew = Xnew_full,
      l_nums = l_nums,
      v_nums = v_nums,
      ...
    )

    is_classification <- !is.null(prediction$class)

    if (is_classification) {
      df <- data.frame(x1 = grid$x1, x2 = grid$x2,
                       class = factor(prediction$class),
                       max_prob = apply(prediction$mean, 1, max))

      if (engine == "ggplot") {
        p1 <- my_2D_plot_fun_class_ggplot("class", "Predicted Classes", df, X_sub, Y, dims = dims)
        p2 <- my_2D_plot_fun_class_ggplot("max_prob", "Maximum Predicted Probability", df, X_sub, Y, classification = FALSE, dims = dims)
      } else {
        p1 <- my_2D_plot_fun_class("class", "Predicted Classes", df, X_sub, Y, dims = dims)
        p2 <- my_2D_plot_fun_class("max_prob", "Maximum Predicted Probability", df, X_sub, Y, classification = FALSE, dims = dims)
      }
      grid.arrange(p1, p2, ncol = 2)
    } else {
      for (j in 1:q) {
        df <- data.frame(x1 = grid$x1, x2 = grid$x2,
                         Mean = prediction$mean[, j],
                         Upper = prediction$upper[, j],
                         Lower = prediction$lower[, j],
                         Variance = prediction$variance[, j])

        if (only_mean) {
          p1 <- if (engine == "ggplot") {
            my_2D_plot_fun_ggplot("Mean", "Predictive Mean", df, dims = dims)
          } else {
            my_2D_plot_fun("Mean", "Predictive Mean", df, dims = dims)
          }
          print(p1)
        } else {
          if (engine == "ggplot") {
            p1 <- my_2D_plot_fun_ggplot("Mean", "Predictive Mean", df, dims = dims)
            p2 <- my_2D_plot_fun_ggplot("Upper", paste0(prediction$CI_level * 100, "% CI Upper"), df, dims = dims)
            p3 <- my_2D_plot_fun_ggplot("Variance", "Predictive Variance", df, dims = dims)
            p4 <- my_2D_plot_fun_ggplot("Lower", paste0(prediction$CI_level * 100, "% CI Lower"), df, dims = dims)
          } else {
            p1 <- my_2D_plot_fun("Mean", "Predictive Mean", df, dims = dims)
            p2 <- my_2D_plot_fun("Upper", paste0(prediction$CI_level * 100, "% CI Upper"), df, dims = dims)
            p3 <- my_2D_plot_fun("Variance", "Predictive Variance", df, dims = dims)
            p4 <- my_2D_plot_fun("Lower", paste0(prediction$CI_level * 100, "% CI Lower"), df, dims = dims)
          }
          grid.arrange(p1, p2, p3, p4, ncol = 2,
                       top = textGrob(paste0("Estimated Probability (Class ", j, ")"),
                                      gp = gpar(fontface = "bold", fontsize = 16)))
        }
      }
    }
  }
}
