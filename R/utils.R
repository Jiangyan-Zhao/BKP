
#' Internal Twin index selector
#' @noRd
get_twin_indices <- function(Xnorm, g, v = 2L * g, runs = 10L, seed = 123L) {
  get_twin_indices_rcpp(
    data = as.matrix(Xnorm),
    g = as.integer(g),
    v = as.integer(v),
    runs = as.integer(runs),
    seed = as.integer(seed)
  )
}

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

  p <- ggplot(data, aes(x = .data$x1, y = .data$x2)) +
    geom_raster(aes(fill = .data[[var]])) +
    geom_contour(aes(z = .data[[var]]), color = "black", linewidth = 0.2) +
    # Remove padding between the plot area and the panel border to mimic base R tightly fitted box
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_viridis_c(
      option = "plasma",
      name = NULL, # Remove the legend title to match base R layout
      limits = fill_limits,
      # Configure colorbar to strictly match base R styling
      guide = guide_colorbar(
        frame.colour = "black",
        frame.linewidth = 0.5, # Sync legend frame thickness with panel border
        ticks.colour = "black",
        # Calculate precise colorbar height:
        # '1 npc' maps to the total plot height. We subtract roughly 5 lines of text
        # (title + axis labels + axis titles + margins) to make the bar exactly match the panel height.
        barheight = unit(0.2, "npc"),
        barwidth = unit(1.2, "lines")
      )
    ) +
    labs(
      title = title,
      x = if (is.null(dims)) "x1" else paste0("x", dims[1]),
      y = if (is.null(dims)) "x2" else paste0("x", dims[2])
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid = element_blank(),
      # Match the plot border thickness directly with the legend frame
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(color = "black"),
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
    p <- p + geom_point(
      data = obs_df,
      aes(x = .data$x1, y = .data$x2, shape = .data$cls),
      color = "red", size = 2, inherit.aes = FALSE
    ) + scale_shape_manual(values = c("0" = 4, "1" = 16), guide = "none")
  }

  p
}

my_2D_plot_fun_class_ggplot <- function(var, title, data, X, Y,
                                        classification = TRUE, dims = NULL, ...) {
  # Validate that var is a single character string
  if (!is.character(var) || length(var) != 1) {
    stop("`var` must be a single character string (a column name).")
  }
  # Ensure the variable exists in the dataset if not in classification mode
  if (!classification && !var %in% names(data)) {
    stop(sprintf("Column `%s` not found in `data`.", var))
  }

  # Identify the class with maximum probability for each observation
  class_Y <- max.col(Y)

  p <- ggplot(data, aes(x = .data$x1, y = .data$x2))

  if (classification) {
    # Categorical background for predicted classes
    p <- p + geom_raster(aes(fill = .data$class), alpha = 0.8) +
      scale_fill_brewer(palette = "Set2", name = NULL)
  } else {
    # Continuous background for probability surfaces (e.g., Max Prob)
    p <- p + geom_raster(aes(fill = .data[[var]])) +
      geom_contour(aes(z = .data[[var]]), color = "black", linewidth = 0.2) +
      scale_fill_viridis_c(
        option = "plasma",
        direction = -1,
        name = NULL,
        guide = guide_colorbar(
          frame.colour = "black",
          frame.linewidth = 0.5,
          ticks.colour = "black",
          barheight = unit(0.5, "npc") # Use relative height to match the panel
        )
      )
  }

  # Prepare observation data for overlay
  obs_df <- data.frame(
    x1 = X[, 1],
    x2 = X[, 2],
    class = factor(class_Y)
  )

  # Overlay observed points with distinct shapes
  p <- p + geom_point(
    data = obs_df,
    aes(x = .data$x1, y = .data$x2, shape = .data$class),
    color = "black", fill = "white", size = 2, stroke = 1,
    inherit.aes = FALSE
  ) +
    # Ensure axes fit tightly to the data grid as in base R
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_shape_manual(values = seq(15, 15 + length(unique(class_Y))), guide = "none") +
    labs(
      title = title,
      x = if (is.null(dims)) "x1" else paste0("x", dims[1]),
      y = if (is.null(dims)) "x2" else paste0("x", dims[2])
    ) +
    # Use standard white background with black border
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.title = element_text(size = 12),
      axis.text  = element_text(color = "black"),
      legend.background = element_blank(),
      legend.key = element_blank()
    )

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


# Beta-binomial helpers (internal)
#' @noRd
betabinom_pmf <- function(k, size, alpha, beta) {
  exp(lchoose(size, k) + lbeta(k + alpha, size - k + beta) - lbeta(alpha, beta))
}

#' @noRd
betabinom_quantile <- function(prob, size, alpha, beta) {
  if (!is.finite(prob) || prob <= 0) return(0)
  if (prob >= 1) return(size)
  k <- 0:size
  cdf <- cumsum(betabinom_pmf(k, size, alpha, beta))
  which(cdf >= prob)[1] - 1L
}
