#' @rdname plot
#'
#' @keywords TwinBKP
#'
#' @export
#' @method plot TwinBKP
plot.TwinBKP <- function(x, only_mean = FALSE, n_grid = 80, dims = NULL,
                         engine = c("base", "ggplot"),
                         show_global = TRUE, ...) {

  # ---------------- Argument Checking ----------------
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
  y <- as.numeric(x$y)
  m <- as.numeric(x$m)
  Xbounds <- x$Xbounds

  d <- ncol(X)

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

  global_indices <- x$global_indices
  has_global <- isTRUE(show_global) &&
    !is.null(global_indices) &&
    length(global_indices) > 0

  if (has_global) {
    X_global_sub <- X[global_indices, dims, drop = FALSE]
    y_global_prop <- y[global_indices] / m[global_indices]
  } else {
    X_global_sub <- NULL
    y_global_prop <- NULL
  }

  # ---------------- 1D Plot ----------------
  if (length(dims) == 1) {

    Xnew <- matrix(
      seq(Xbounds[dims, 1], Xbounds[dims, 2], length.out = 10 * n_grid),
      ncol = 1
    )

    Xnew_full <- make_plot_grid(X, dims, Xnew)

    prediction <- predict.TwinBKP(x, Xnew = Xnew_full, ...)

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
      lbl_ci <- paste0(prediction$CI_level * 100, "% Credible Interval")
      lbl_pts <- "Observed Proportions"

      p <- ggplot(plot_df, aes(x = .data$x)) +
        geom_ribbon(
          aes(ymin = .data$lower, ymax = .data$upper),
          fill = "grey70",
          alpha = 0.4
        ) +
        geom_line(
          aes(y = .data$mean, color = lbl_ci),
          alpha = 0
        ) +
        geom_line(
          aes(y = .data$mean, color = lbl_line),
          linewidth = 1
        ) +
        geom_point(
          data = obs_df,
          aes(x = .data$x, y = .data$obs_y, color = lbl_pts),
          size = 2
        ) +
        scale_color_manual(
          name = NULL,
          values = stats::setNames(
            c("blue", "grey70", "red"),
            c(lbl_line, lbl_ci, lbl_pts)
          ),
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
          title = "TwinBKP Estimated Probability",
          x = ifelse(d > 1, paste0("x", dims), "x"),
          y = "Probability"
        ) +
        theme_bw() +
        theme(
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 11, color = "black"),
          legend.position = c(0.02, 0.98),
          legend.justification = c(0, 1),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(size = 12),
          legend.key.width = unit(2, "line")
        )

      if (has_global) {
        global_df <- data.frame(
          x = as.numeric(X_global_sub),
          obs_y = y_global_prop
        )

        p <- p +
          geom_point(
            data = global_df,
            aes(x = .data$x, y = .data$obs_y),
            inherit.aes = FALSE,
            shape = 21,
            fill = "white",
            color = "black",
            size = 2.8,
            stroke = 0.9
          )
      }

      if (is_classification) {
        p <- p +
          geom_hline(yintercept = prediction$threshold, linetype = "dashed") +
          annotate(
            "text",
            x = max(plot_df$x),
            y = prediction$threshold + 0.02,
            label = "threshold",
            hjust = 1,
            vjust = 0.5
          )
      }

      print(p)

    } else {

      plot(
        Xnew, prediction$mean,
        type = "l",
        col = "blue",
        lwd = 2,
        xlab = ifelse(d > 1, paste0("x", dims), "x"),
        ylab = "Probability",
        main = "TwinBKP Estimated Probability",
        xlim = Xbounds[dims, ],
        ylim = c(
          max(0, min(prediction$lower) * 0.9),
          min(1, max(prediction$upper) * 1.1)
        )
      )

      polygon(
        c(Xnew, rev(Xnew)),
        c(prediction$lower, rev(prediction$upper)),
        col = "lightgrey",
        border = NA
      )

      lines(Xnew, prediction$mean, col = "blue", lwd = 2)

      points(X_sub, y / m, pch = 20, col = "red")

      if (has_global) {
        points(
          X_global_sub,
          y_global_prop,
          pch = 21,
          bg = "white",
          col = "black",
          cex = 1.25,
          lwd = 1
        )
      }

      if (is_classification) {
        abline(h = prediction$threshold, lty = 2, lwd = 1.2)
        text(
          x = Xbounds[dims, 2],
          y = prediction$threshold + 0.02,
          labels = "threshold",
          adj = c(1, 0.5),
          cex = 0.9,
          col = "black"
        )
      }

      legend_items <- c(
        "Estimated Probability",
        paste0(prediction$CI_level * 100, "% Credible Interval"),
        "Observed Proportions"
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

    return(invisible(NULL))
  }

  # ---------------- 2D Plot ----------------
  x1 <- seq(Xbounds[dims[1], 1], Xbounds[dims[1], 2], length.out = n_grid)
  x2 <- seq(Xbounds[dims[2], 1], Xbounds[dims[2], 2], length.out = n_grid)
  grid <- expand.grid(x1 = x1, x2 = x2)

  Xnew_full <- make_plot_grid(X, dims, grid)

  prediction <- predict.TwinBKP(x, Xnew = Xnew_full, ...)

  is_classification <- !is.null(prediction$class)

  df <- data.frame(
    x1 = grid$x1,
    x2 = grid$x2,
    Mean = prediction$mean,
    Upper = prediction$upper,
    Lower = prediction$lower,
    Variance = prediction$variance
  )

  if (has_global) {
    X_global_plot <- X_global_sub
  } else {
    X_global_plot <- NULL
  }

  make_plot <- function(var, title, overlay_obs = FALSE) {
    if (engine == "ggplot") {
      my_2D_plot_fun_ggplot(
        var = var,
        title = title,
        data = df,
        X = if (overlay_obs) X_sub else NULL,
        y = if (overlay_obs) y else NULL,
        X_global = X_global_plot,
        dims = dims
      )
    } else {
      my_2D_plot_fun(
        var = var,
        title = title,
        data = df,
        X = if (overlay_obs) X_sub else NULL,
        y = if (overlay_obs) y else NULL,
        X_global = X_global_plot,
        dims = dims
      )
    }
  }

  if (only_mean) {
    if (is_classification) {
      p1 <- make_plot(
        "Mean",
        "Predicted Class Probability (Predictive Mean)",
        overlay_obs = TRUE
      )
    } else {
      p1 <- make_plot(
        "Mean",
        "TwinBKP Predictive Mean",
        overlay_obs = FALSE
      )
    }

    print(p1)
    return(invisible(NULL))
  }

  if (is_classification) {
    p1 <- make_plot("Mean", "TwinBKP Predictive Mean", overlay_obs = TRUE)
    p3 <- make_plot("Variance", "TwinBKP Predictive Variance", overlay_obs = TRUE)

    gridExtra::grid.arrange(p1, p3, ncol = 2)
  } else {
    p1 <- make_plot("Mean", "TwinBKP Predictive Mean", overlay_obs = FALSE)
    p2 <- make_plot(
      "Upper",
      paste0(prediction$CI_level * 100, "% CI Upper"),
      overlay_obs = FALSE
    )
    p3 <- make_plot("Variance", "TwinBKP Predictive Variance", overlay_obs = FALSE)
    p4 <- make_plot(
      "Lower",
      paste0(prediction$CI_level * 100, "% CI Lower"),
      overlay_obs = FALSE
    )

    gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
  }

  return(invisible(NULL))
}
