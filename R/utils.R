
# Helper to build each plot
#' @noRd

my_2D_plot_fun <- function(var, title, data, X = NULL, y = NULL, dims = NULL, ...) {
  levelplot(
    as.formula(paste(var, "~ x1 * x2")),
    data = data,
    col.regions = hcl.colors(100, palette = "plasma"),
    main = title,
    xlab = ifelse(is.null(dims), "x1", paste0("x", dims[1])),
    ylab = ifelse(is.null(dims), "x2", paste0("x", dims[2])),
    contour = TRUE,
    colorkey = TRUE,
    cuts = 15,
    pretty = TRUE,
    scales = list(draw = TRUE, tck = c(1, 0)),
    panel = function(...) {
      panel.levelplot(...)
      panel.contourplot(..., col = "black", lwd = 0.5)
      panel.points(X[,1], X[,2], pch = ifelse(y == 1, 16, 4),
                   col = "red", lwd = 2, cex = 1.2)
    }
  )
}


my_2D_plot_fun_class <- function(var, title, data, X, Y, classification = TRUE, dims = NULL, ...) {
  class_Y <- max.col(Y)

  if(classification){
    q <- ncol(Y)
    cols <- hcl.colors(q, palette = "Cold")
    colorkey <- FALSE
    cuts <- q
  }else{
    cols <- hcl.colors(100, palette = "plasma", rev = TRUE)
    colorkey <- TRUE
    cuts <- 15
  }

  levelplot(
    as.formula(paste(var, "~ x1 * x2")),
    data = data,
    col.regions = cols,
    main = title,
    xlab = ifelse(is.null(dims), "x1", paste0("x", dims[1])),
    ylab = ifelse(is.null(dims), "x2", paste0("x", dims[2])),
    colorkey = colorkey,
    cuts = cuts,
    pretty = TRUE,
    scales = list(draw = TRUE, tck = c(1, 0)),
    panel = function(...) {
      panel.levelplot(...)
      panel.contourplot(..., col = "black", lwd = 0.5)
      panel.points(X[, 1], X[, 2], pch = class_Y, col = "black",
                   fill = cols[class_Y], lwd = 1.5, cex = 1.2)
    }
  )
}

my_2D_plot_fun_ggplot <- function(var, title, data, X = NULL, y = NULL, dims = NULL, ...) {
  if (!is.character(var) || length(var) != 1) {
    stop("`var` must be a single character string (a column name).")
  }
  if (!var %in% names(data)) {
    stop(sprintf("Column `%s` not found in `data`.", var))
  }

  p <- ggplot(data, aes(x = .data$x1, y = .data$x2)) +
    geom_raster(aes(fill = .data[[var]])) +
    geom_contour(aes(z = .data[[var]]), color = "black", linewidth = 0.2) +
    scale_fill_viridis_c(option = "plasma") +
    labs(
      title = title,
      x = if (is.null(dims)) "x1" else paste0("x", dims[1]),
      y = if (is.null(dims)) "x2" else paste0("x", dims[2]),
      fill = var
    ) + theme_minimal()

  if (!is.null(X) && !is.null(y)) {
    obs_df <- data.frame(
      x1 = X[, 1],
      x2 = X[, 2],
      cls = factor(ifelse(y == 1, "1", "0"))
    )
    p <- p + geom_point(
        data = obs_df,
        aes(x = .data$x1, y = .data$x2, shape = .data$cls),
        color = "red", size = 2, inherit.aes = FALSE
      ) + scale_shape_manual(values = c("0" = 4, "1" = 16), guide = "none")
  }

  p
}

my_2D_plot_fun_ggplot <- function(var, title, data, X = NULL, y = NULL, dims = NULL, ...) {
  # Validate inputs to ensure 'var' is a valid column string
  if (!is.character(var) || length(var) != 1) {
    stop("`var` must be a single character string (a column name).")
  }
  if (!var %in% names(data)) {
    stop(sprintf("Column `%s` not found in `data`.", var))
  }

  # Constrain probability metrics (Mean, Upper, Lower) to a strict [0, 1] range.
  # Variance is left unconstrained (NULL) to preserve color gradients for small values.
  fill_limits <- if (var %in% c("Mean", "Upper", "Lower")) c(0, 1) else NULL

  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$x1, y = .data$x2)) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data[[var]])) +
    ggplot2::geom_contour(ggplot2::aes(z = .data[[var]]), color = "black", linewidth = 0.2) +
    # Remove padding between the plot area and the panel border to mimic base R tightly fitted box
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_fill_viridis_c(
      option = "plasma",
      name = NULL, # Remove the legend title to match base R layout
      limits = fill_limits,
      # Configure colorbar to strictly match base R styling
      guide = ggplot2::guide_colorbar(
        frame.colour = "black",
        frame.linewidth = 0.5, # Sync legend frame thickness with panel border
        ticks.colour = "black",
        # Calculate precise colorbar height:
        # '1 npc' maps to the total plot height. We subtract roughly 5 lines of text
        # (title + axis labels + axis titles + margins) to make the bar exactly match the panel height.
        barheight = ggplot2::unit(0.2, "npc"),
        barwidth = ggplot2::unit(1.2, "lines")
      )
    ) +
    ggplot2::labs(
      title = title,
      x = if (is.null(dims)) "x1" else paste0("x", dims[1]),
      y = if (is.null(dims)) "x2" else paste0("x", dims[2])
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid = ggplot2::element_blank(),
      # Match the plot border thickness directly with the legend frame
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(color = "black"),
      # Center the legend vertically to align perfectly with the panel
      legend.justification = "center"
    )

  # Overlay observation points conditionally
  if (!is.null(X) && !is.null(y)) {
    obs_df <- data.frame(
      x1 = X[, 1],
      x2 = X[, 2],
      cls = factor(ifelse(y == 1, "1", "0"))
    )
    p <- p + ggplot2::geom_point(
      data = obs_df,
      ggplot2::aes(x = .data$x1, y = .data$x2, shape = .data$cls),
      color = "red", size = 2, inherit.aes = FALSE
    ) + ggplot2::scale_shape_manual(values = c("0" = 4, "1" = 16), guide = "none")
  }

  p
}
posterior_summary <- function(mean_vals, var_vals) {
  summary_mat <- rbind(
    "Posterior means" = c(
      Mean   = mean(mean_vals),
      Median = median(mean_vals),
      SD     = sd(mean_vals),
      Min    = min(mean_vals),
      Max    = max(mean_vals)
    ),
    "Posterior variances" = c(
      Mean   = mean(var_vals),
      Median = median(var_vals),
      SD     = sd(var_vals),
      Min    = min(var_vals),
      Max    = max(var_vals)
    )
  )
  return(round(summary_mat, 4))
}
