
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
