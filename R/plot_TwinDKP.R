#' @rdname plot
#'
#' @keywords TwinDKP
#'
#' @export
#' @method plot TwinDKP
plot.TwinDKP <- function(
    x,
    only_mean = FALSE,
    n_grid = 80,
    dims = NULL,
    engine = c("base", "ggplot"),
    show_global = TRUE,
    ...
) {
  # ---------------- Argument Checking ----------------
  if (!inherits(x, "TwinDKP")) {
    stop("'x' must be a fitted TwinDKP object.")
  }

  if (!is.logical(only_mean) || length(only_mean) != 1) {
    stop("`only_mean` must be a single logical value (TRUE or FALSE).")
  }

  if (!is.numeric(n_grid) || length(n_grid) != 1 ||
      n_grid <= 0 || n_grid != floor(n_grid)) {
    stop("'n_grid' must be a positive integer.")
  }
  n_grid <- as.integer(n_grid)

  if (!is.logical(show_global) || length(show_global) != 1) {
    stop("'show_global' must be a single logical value.")
  }

  engine <- match.arg(engine)

  # ---------------- Extract Model Components ----------------
  X <- x$X
  Y <- x$Y
  Xbounds <- x$Xbounds

  d <- ncol(X)
  q <- ncol(Y)

  # ---------------- Handle dims Argument ----------------
  if (is.null(dims)) {
    if (d > 2) {
      stop("X has more than 2 dimensions. Please specify `dims` for plotting.")
    }
    dims <- seq_len(d)
  } else {
    if (!is.numeric(dims) || any(dims != as.integer(dims))) {
      stop("`dims` must be an integer vector.")
    }
    dims <- as.integer(dims)

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

  X_sub <- X[, dims, drop = FALSE]
  is_classification <- all(rowSums(Y) == 1)
  class_Y <- if (is_classification) max.col(Y) else NULL

  global_indices <- x$global_indices
  has_global <- isTRUE(show_global) &&
    !is.null(global_indices) &&
    length(global_indices) > 0

  if (has_global) {
    X_global_sub <- X[global_indices, dims, drop = FALSE]
    Y_global <- Y[global_indices, , drop = FALSE]
    Y_global_prop <- sweep(Y_global, 1, rowSums(Y_global), "/")
    class_global <- if (is_classification) max.col(Y_global) else NULL
  } else {
    X_global_sub <- NULL
    Y_global_prop <- NULL
    class_global <- NULL
  }

  # ---------------- 1D Plot ----------------
  if (length(dims) == 1) {
    Xnew <- matrix(
      seq(Xbounds[dims, 1], Xbounds[dims, 2], length.out = 10 * n_grid),
      ncol = 1
    )

    Xnew_full <- .make_plot_grid(X, Xbounds, dims, Xnew)
    prediction <- predict.TwinDKP(x, Xnew = Xnew_full, ...)

    if (engine == "ggplot") {
      plot_list <- vector("list", q + is_classification)

      lbl_line <- "Estimated Probability"
      lbl_ci <- paste0(prediction$CI_level * 100, "% CI")
      lbl_pts <- "Observed"

      if (is_classification) {
        all_pred_df <- data.frame(
          x = rep(as.numeric(Xnew), q),
          prob = as.vector(prediction$mean),
          Class = factor(rep(seq_len(q), each = nrow(Xnew)))
        )
        all_obs_df <- data.frame(
          x = as.numeric(X_sub),
          y = rep(-0.05, nrow(X_sub)),
          Class = factor(class_Y, levels = seq_len(q))
        )

        p_all <- ggplot() +
          geom_line(
            data = all_pred_df,
            aes(x = .data$x, y = .data$prob, color = .data$Class),
            linewidth = 1
          ) +
          geom_point(
            data = all_obs_df,
            aes(x = .data$x, y = .data$y, color = .data$Class),
            size = 2
          ) +
          scale_color_discrete(name = NULL, labels = paste("Class", seq_len(q))) +
          labs(
            title = "TwinDKP Estimated Mean Curves (All Classes)",
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

        if (has_global) {
          global_all_df <- data.frame(
            x = as.numeric(X_global_sub),
            y = rep(-0.05, nrow(X_global_sub)),
            Class = factor(class_global, levels = seq_len(q))
          )
          p_all <- p_all +
            geom_point(
              data = global_all_df,
              aes(x = .data$x, y = .data$y),
              inherit.aes = FALSE,
              shape = 21,
              fill = "white",
              color = "black",
              size = 2.8,
              stroke = 0.9
            )
        }

        plot_list[[1]] <- p_all
      }

      for (j in seq_len(q)) {
        mean_j <- prediction$mean[, j]
        lower_j <- prediction$lower[, j]
        upper_j <- prediction$upper[, j]

        pred_df_j <- data.frame(
          x = as.numeric(Xnew),
          mean = mean_j,
          lower = lower_j,
          upper = upper_j
        )

        if (is_classification) {
          obs_j <- as.integer(class_Y == j)
          global_obs_j <- as.integer(class_global == j)
          ylim_j <- c(0, 1)
        } else {
          obs_j <- Y[, j] / rowSums(Y)
          global_obs_j <- if (has_global) Y_global_prop[, j] else NULL
          ylim_j <- c(max(0, min(lower_j) * 0.9), min(1, max(upper_j) * 1.1))
        }

        obs_df_j <- data.frame(x = as.numeric(X_sub), obs = obs_j)

        p <- ggplot() +
          geom_ribbon(
            data = pred_df_j,
            aes(x = .data$x, ymin = .data$lower, ymax = .data$upper),
            fill = "grey70",
            alpha = 0.4
          ) +
          geom_line(data = pred_df_j, aes(x = .data$x, y = .data$mean, color = lbl_ci), alpha = 0) +
          geom_line(data = pred_df_j, aes(x = .data$x, y = .data$mean, color = lbl_line), linewidth = 1) +
          geom_point(data = obs_df_j, aes(x = .data$x, y = .data$obs, color = lbl_pts), size = 2) +
          scale_color_manual(
            name = NULL,
            values = stats::setNames(c("blue", "grey70", "red"), c(lbl_line, lbl_ci, lbl_pts)),
            breaks = c(lbl_line, lbl_ci, lbl_pts)
          ) +
          guides(
            color = guide_legend(
              override.aes = list(
                shape = c(NA, NA, 16),
                linetype = c(1, 1, 0),
                linewidth = c(1, 5, 0),
                alpha = c(1, 0.5, 1)
              )
            )
          ) +
          labs(
            title = paste0("TwinDKP Estimated Probability (Class ", j, ")"),
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
            axis.text = element_text(size = 10, color = "black")
          )

        if (has_global) {
          global_df_j <- data.frame(
            x = as.numeric(X_global_sub),
            obs = global_obs_j
          )
          p <- p +
            geom_point(
              data = global_df_j,
              aes(x = .data$x, y = .data$obs),
              inherit.aes = FALSE,
              shape = 21,
              fill = "white",
              color = "black",
              size = 2.8,
              stroke = 0.9
            )
        }

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
      on.exit(par(old_par), add = TRUE)

      if (is_classification) {
        cols <- rainbow(q)
        plot(
          NA,
          xlim = Xbounds[dims, ],
          ylim = c(-0.1, 1.1),
          xlab = ifelse(d > 1, paste0("x", dims), "x"),
          ylab = "Probability",
          main = "TwinDKP Estimated Mean Curves (All Classes)"
        )
        for (j in seq_len(q)) {
          lines(Xnew, prediction$mean[, j], col = cols[j], lwd = 2)
        }
        points(X_sub, rep(-0.05, nrow(X_sub)), col = cols[class_Y], pch = 20)
        if (has_global) {
          points(
            X_global_sub,
            rep(-0.05, nrow(X_global_sub)),
            pch = 21,
            bg = "white",
            col = "black",
            cex = 1.25,
            lwd = 1
          )
        }
        legend(
          "top",
          legend = paste("Class", seq_len(q)),
          col = cols,
          lty = 1,
          lwd = 2,
          horiz = TRUE,
          bty = "n"
        )
      }

      for (j in seq_len(q)) {
        mean_j <- prediction$mean[, j]
        lower_j <- prediction$lower[, j]
        upper_j <- prediction$upper[, j]

        if (is_classification) {
          ylim <- c(0, 1)
          obs_j <- as.integer(class_Y == j)
          global_obs_j <- as.integer(class_global == j)
        } else {
          ylim <- c(max(0, min(lower_j) * 0.9), min(1, max(upper_j) * 1.1))
          obs_j <- Y[, j] / rowSums(Y)
          global_obs_j <- if (has_global) Y_global_prop[, j] else NULL
        }

        plot(
          Xnew,
          mean_j,
          type = "l",
          col = "blue",
          lwd = 2,
          xlab = ifelse(d > 1, paste0("x", dims), "x"),
          ylab = "Probability",
          main = paste0("TwinDKP Estimated Probability (Class ", j, ")"),
          xlim = Xbounds[dims, ],
          ylim = ylim
        )

        polygon(
          c(Xnew, rev(Xnew)),
          c(lower_j, rev(upper_j)),
          col = "lightgrey",
          border = NA
        )
        lines(Xnew, mean_j, col = "blue", lwd = 2)
        points(X_sub, obs_j, pch = 20, col = "red")

        if (has_global) {
          points(
            X_global_sub,
            global_obs_j,
            pch = 21,
            bg = "white",
            col = "black",
            cex = 1.25,
            lwd = 1
          )
        }

        if (j == 1) {
          legend_items <- c(
            "Estimated Probability",
            paste0(prediction$CI_level * 100, "% CI"),
            "Observed"
          )
          legend_col <- c("blue", "lightgrey", "red")
          legend_lwd <- c(2, 8, NA)
          legend_pch <- c(NA, NA, 20)
          legend_pt_bg <- c(NA, NA, NA)

          if (has_global) {
            legend_items <- c(legend_items, "Global Subset")
            legend_col <- c(legend_col, "black")
            legend_lwd <- c(legend_lwd, NA)
            legend_pch <- c(legend_pch, 21)
            legend_pt_bg <- c(legend_pt_bg, "white")
          }

          legend(
            "topleft",
            legend = legend_items,
            col = legend_col,
            lwd = legend_lwd,
            pch = legend_pch,
            pt.bg = legend_pt_bg,
            lty = c(1, 1, NA, if (has_global) NA),
            bty = "n"
          )
        }
      }
    }

    return(invisible(NULL))
  }

  # ---------------- 2D Plot ----------------
  x1 <- seq(Xbounds[dims[1], 1], Xbounds[dims[1], 2], length.out = n_grid)
  x2 <- seq(Xbounds[dims[2], 1], Xbounds[dims[2], 2], length.out = n_grid)
  grid <- expand.grid(x1 = x1, x2 = x2)

  Xnew_full <- .make_plot_grid(X, Xbounds, dims, grid)
  prediction <- predict.TwinDKP(x, Xnew = Xnew_full, ...)

  if (is_classification) {
    class_pred <- max.col(prediction$mean)
    max_prob <- apply(prediction$mean, 1, max)
    df <- data.frame(
      x1 = grid$x1,
      x2 = grid$x2,
      class = factor(class_pred),
      max_prob = max_prob
    )

    if (engine == "ggplot") {
      p1 <- my_2D_plot_fun_class_ggplot(
        "class",
        "Predicted Classes",
        df,
        X_sub,
        Y,
        dims = dims,
        X_global = if (has_global) X_global_sub else NULL
      )
      p2 <- my_2D_plot_fun_class_ggplot(
        "max_prob",
        "Maximum Predicted Probability",
        df,
        X_sub,
        Y,
        classification = FALSE,
        dims = dims,
        X_global = if (has_global) X_global_sub else NULL
      )
    } else {
      p1 <- my_2D_plot_fun_class(
        "class",
        "Predicted Classes",
        df,
        X_sub,
        Y,
        dims = dims,
        X_global = if (has_global) X_global_sub else NULL
      )
      p2 <- my_2D_plot_fun_class(
        "max_prob",
        "Maximum Predicted Probability",
        df,
        X_sub,
        Y,
        classification = FALSE,
        dims = dims,
        X_global = if (has_global) X_global_sub else NULL
      )
    }
    grid.arrange(p1, p2, ncol = 2)
  } else {
    for (j in seq_len(q)) {
      df <- data.frame(
        x1 = grid$x1,
        x2 = grid$x2,
        Mean = prediction$mean[, j],
        Upper = prediction$upper[, j],
        Lower = prediction$lower[, j],
        Variance = prediction$variance[, j]
      )

      if (only_mean) {
        p1 <- if (engine == "ggplot") {
          my_2D_plot_fun_ggplot(
            "Mean",
            "TwinDKP Predictive Mean",
            df,
            dims = dims,
            X_global = if (has_global) X_global_sub else NULL
          )
        } else {
          my_2D_plot_fun(
            "Mean",
            "TwinDKP Predictive Mean",
            df,
            dims = dims,
            X_global = if (has_global) X_global_sub else NULL
          )
        }
        print(p1)
      } else {
        if (engine == "ggplot") {
          p1 <- my_2D_plot_fun_ggplot(
            "Mean",
            "TwinDKP Predictive Mean",
            df,
            dims = dims,
            X_global = if (has_global) X_global_sub else NULL
          )
          p2 <- my_2D_plot_fun_ggplot(
            "Upper",
            paste0(prediction$CI_level * 100, "% CI Upper"),
            df,
            dims = dims,
            X_global = if (has_global) X_global_sub else NULL
          )
          p3 <- my_2D_plot_fun_ggplot(
            "Variance",
            "TwinDKP Predictive Variance",
            df,
            dims = dims,
            X_global = if (has_global) X_global_sub else NULL
          )
          p4 <- my_2D_plot_fun_ggplot(
            "Lower",
            paste0(prediction$CI_level * 100, "% CI Lower"),
            df,
            dims = dims,
            X_global = if (has_global) X_global_sub else NULL
          )
        } else {
          p1 <- my_2D_plot_fun(
            "Mean",
            "TwinDKP Predictive Mean",
            df,
            dims = dims,
            X_global = if (has_global) X_global_sub else NULL
          )
          p2 <- my_2D_plot_fun(
            "Upper",
            paste0(prediction$CI_level * 100, "% CI Upper"),
            df,
            dims = dims,
            X_global = if (has_global) X_global_sub else NULL
          )
          p3 <- my_2D_plot_fun(
            "Variance",
            "TwinDKP Predictive Variance",
            df,
            dims = dims,
            X_global = if (has_global) X_global_sub else NULL
          )
          p4 <- my_2D_plot_fun(
            "Lower",
            paste0(prediction$CI_level * 100, "% CI Lower"),
            df,
            dims = dims,
            X_global = if (has_global) X_global_sub else NULL
          )
        }

        grid.arrange(
          p1,
          p2,
          p3,
          p4,
          ncol = 2,
          top = textGrob(
            paste0("TwinDKP Estimated Probability (Class ", j, ")"),
            gp = gpar(fontface = "bold", fontsize = 16)
          )
        )
      }
    }
  }

  invisible(NULL)
}
