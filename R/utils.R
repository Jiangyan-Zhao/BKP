
# Helper to build each plot
#' @noRd

my_2D_plot_fun <- function(var, title, data, X = NULL, y = NULL, ...) {
  stopifnot(ncol(X) == 2)
  levelplot(
    as.formula(paste(var, "~ x1 * x2")),
    data = data,
    col.regions = hcl.colors(100, palette = "plasma"),
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


my_2D_plot_fun_class <- function(var, title, data, X, Y, classification = TRUE, ...) {
  stopifnot(ncol(X) == 2)

  class_Y <- max.col(Y)

  if(classification){
    q <- ncol(Y)
    cols <- hcl.colors(q, palette = "Set2")
    cuts <- q
  }else{
    cols <- hcl.colors(100, palette = "plasma", rev = TRUE)
    cuts <- 15
  }

  levelplot(
    as.formula(paste(var, "~ x1 * x2")),
    data = data,
    col.regions = cols,
    main = title,
    xlab = "X1", ylab = "X2",
    colorkey = TRUE,
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
