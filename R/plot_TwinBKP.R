#' @rdname plot
#' @keywords BKP
#'
#' @details
#' \code{plot.TwinBKP()} follows the same high-level workflow as
#' \code{plot.BKP()}, but all predictive quantities are produced through
#' \code{predict.TwinBKP()} on a constructed grid.
#'
#' For 1D: a dense line grid is generated, then mean and credible interval are
#' plotted with observed proportions.
#'
#' For 2D: a lattice/ggplot surface grid is generated, then mean/upper/lower/
#' variance panels are visualized from prediction outputs.
#'
#' @param l_nums Optional local neighbor count passed to
#'   \code{predict.TwinBKP()}.
#' @param v_nums Optional validation size passed to
#'   \code{predict.TwinBKP()}.
#'
#' @export
#' @method plot TwinBKP
plot.TwinBKP <- function(x, only_mean = FALSE, n_grid = 80, dims = NULL,
                         engine = c("base", "ggplot"),
                         l_nums = NULL, v_nums = NULL, ...) {
  if (!is.logical(only_mean) || length(only_mean) != 1) {
    stop("`only_mean` must be a single logical value (TRUE or FALSE).")
  }
  if (!is.numeric(n_grid) || length(n_grid) != 1 || n_grid <= 0) {
    stop("'n_grid' must be a positive integer.")
  }
  n_grid <- as.integer(n_grid)
  engine <- match.arg(engine)

  X <- x$X
  y <- x$y
  m <- x$m
  Xbounds <- x$Xbounds
  d <- ncol(X)

  if (is.null(dims)) {
    if (d > 2) stop("X has more than 2 dimensions. Please specify `dims` for plotting.")
    dims <- seq_len(d)
  }
  if (length(dims) < 1 || length(dims) > 2) stop("`dims` must have length 1 or 2.")

  X_sub <- X[, dims, drop = FALSE]

  if (length(dims) == 1) {
    Xnew <- matrix(seq(Xbounds[dims, 1], Xbounds[dims, 2], length.out = 10 * n_grid), ncol = 1)
    Xnew_full <- tgp::lhs(nrow(Xnew), Xbounds)
    Xnew_full[, dims] <- Xnew
    pred <- predict.TwinBKP(x, Xnew = Xnew_full, l_nums = l_nums, v_nums = v_nums, ...)

    if (engine == "ggplot") {
      plot_df <- data.frame(
        x = as.numeric(Xnew),
        mean = as.numeric(pred$mean),
        lower = as.numeric(pred$lower),
        upper = as.numeric(pred$upper)
      )
      obs_df <- data.frame(x = as.numeric(X_sub), obs_y = as.numeric(y / m))

      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                             fill = "grey70", alpha = 0.4) +
        ggplot2::geom_line(ggplot2::aes(y = .data$mean), color = "blue", linewidth = 1) +
        ggplot2::geom_point(data = obs_df, ggplot2::aes(y = .data$obs_y), color = "red", size = 2) +
        ggplot2::labs(title = "TwinBKP Estimated Probability",
                      x = ifelse(d > 1, paste0("x", dims), "x"), y = "Probability") +
        ggplot2::theme_bw()
      print(p)
    } else {
      plot(Xnew, as.numeric(pred$mean),
           type = "l", col = "blue", lwd = 2,
           xlab = ifelse(d > 1, paste0("x", dims), "x"),
           ylab = "Probability",
           main = "TwinBKP Estimated Probability",
           xlim = Xbounds[dims, ], ylim = c(0, 1))
      polygon(c(Xnew, rev(Xnew)),
              c(as.numeric(pred$lower), rev(as.numeric(pred$upper))),
              col = "lightgrey", border = NA)
      lines(Xnew, as.numeric(pred$mean), col = "blue", lwd = 2)
      points(X_sub, y / m, pch = 20, col = "red")
    }
  } else {
    x1 <- seq(Xbounds[dims[1], 1], Xbounds[dims[1], 2], length.out = n_grid)
    x2 <- seq(Xbounds[dims[2], 1], Xbounds[dims[2], 2], length.out = n_grid)
    grid <- expand.grid(x1 = x1, x2 = x2)

    Xnew_full <- tgp::lhs(nrow(grid), Xbounds)
    Xnew_full[, dims] <- as.matrix(grid)
    pred <- predict.TwinBKP(x, Xnew = Xnew_full, l_nums = l_nums, v_nums = v_nums, ...)

    df <- data.frame(
      x1 = grid$x1,
      x2 = grid$x2,
      Mean = as.numeric(pred$mean),
      Upper = as.numeric(pred$upper),
      Lower = as.numeric(pred$lower),
      Variance = as.numeric(pred$variance)
    )

    if (only_mean) {
      p1 <- if (engine == "ggplot") {
        my_2D_plot_fun_ggplot("Mean", "TwinBKP Predictive Mean", df, dims = dims)
      } else {
        my_2D_plot_fun("Mean", "TwinBKP Predictive Mean", df, dims = dims)
      }
      print(p1)
    } else {
      if (engine == "ggplot") {
        p1 <- my_2D_plot_fun_ggplot("Mean", "TwinBKP Predictive Mean", df, dims = dims)
        p2 <- my_2D_plot_fun_ggplot("Upper", paste0(pred$CI_level * 100, "% CI Upper"), df, dims = dims)
        p3 <- my_2D_plot_fun_ggplot("Variance", "TwinBKP Predictive Variance", df, dims = dims)
        p4 <- my_2D_plot_fun_ggplot("Lower", paste0(pred$CI_level * 100, "% CI Lower"), df, dims = dims)
      } else {
        p1 <- my_2D_plot_fun("Mean", "TwinBKP Predictive Mean", df, dims = dims)
        p2 <- my_2D_plot_fun("Upper", paste0(pred$CI_level * 100, "% CI Upper"), df, dims = dims)
        p3 <- my_2D_plot_fun("Variance", "TwinBKP Predictive Variance", df, dims = dims)
        p4 <- my_2D_plot_fun("Lower", paste0(pred$CI_level * 100, "% CI Lower"), df, dims = dims)
      }
      gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
    }
  }
}
