
# Helper to build each plot
#' @noRd

my_2D_plot_fun <- function(var, title, data, pal = "plasma", X = NULL, y = NULL, ...) {
  levelplot(
    as.formula(paste(var, "~ x1 * x2")),
    data = data,
    col.regions = hcl.colors(100, palette = pal),
    main = title,
    xlab = "X1", ylab = "X2",
    contour = TRUE,
    colorkey = TRUE,
    cuts = 15,
    pretty = TRUE,
    scales = list(draw = TRUE, tck = c(1, 0)),
    panel = function(...) {
      panel.levelplot(...)
      panel.contourplot(..., col = "black", lwd = 0.5)
      panel.points(X[,1], X[,2], pch = ifelse(y == 1, 16, 15), cex = 1.2)
    }
  )
}
