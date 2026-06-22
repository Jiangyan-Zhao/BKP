#' Create a 2D lattice contour plot for BKP posterior summaries
#'
#' Internal helper used by \code{plot.BKP()} and related plotting methods to
#' create lattice-based contour plots for two-dimensional posterior summaries,
#' such as the predictive mean, variance, and credible interval bounds. Optional
#' binary observations can be overlaid on the surface.
#'
#' @param var Character string giving the name of the column in \code{data} to
#'   plot.
#' @param title Character string used as the plot title.
#' @param data Data frame containing columns \code{x1}, \code{x2}, and the
#'   summary variable named by \code{var}.
#' @param X Optional numeric matrix of observed two-dimensional input locations
#'   to overlay.
#' @param y Optional binary response vector corresponding to \code{X}. Points
#'   with \code{y == 1} and \code{y != 1} are drawn with different symbols.
#' @param dims Optional integer vector of length two indicating the original
#'   input dimensions being visualized. Used only for axis labels.
#' @param ... Additional arguments reserved for future use.
#'
#' @return A \code{trellis} object produced by \code{lattice::levelplot()}.
#'
#' @keywords internal

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
      if (!is.null(X) && !is.null(y)) {
        lattice::panel.points(
          X[, 1], X[, 2],
          pch = ifelse(y == 1, 16, 4),
          col = "red", lwd = 2, cex = 1.2
        )
      }
    }
  )
}


#' Create a 2D lattice plot for DKP class predictions or probabilities
#'
#' Internal helper used by \code{plot.DKP()} to create lattice-based two-
#' dimensional plots for predicted classes or maximum posterior class
#' probabilities. Observed class labels are overlaid using symbols determined by
#' the observed multinomial counts.
#'
#' @param var Character string giving the name of the column in \code{data} to
#'   plot.
#' @param title Character string used as the plot title.
#' @param data Data frame containing columns \code{x1}, \code{x2}, and the
#'   variable named by \code{var}.
#' @param X Numeric matrix of observed two-dimensional input locations to
#'   overlay.
#' @param Y Numeric matrix of multinomial counts or one-hot class indicators.
#'   The observed class is determined by \code{max.col(Y)}.
#' @param classification Logical. If \code{TRUE}, the plotted surface is treated
#'   as categorical class output; otherwise it is treated as a continuous
#'   probability or uncertainty surface.
#' @param dims Optional integer vector of length two indicating the original
#'   input dimensions being visualized. Used only for axis labels.
#' @param ... Additional arguments reserved for future use.
#'
#' @return A \code{trellis} object produced by \code{lattice::levelplot()}.
#'
#' @keywords internal

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


#' Create a styled 2D ggplot surface for BKP posterior summaries
#'
#' Internal helper used by \code{plot.BKP()} to create a styled
#' \pkg{ggplot2}-based raster and contour plot for two-dimensional posterior
#' summaries.
#'
#' @param var Character string giving the name of the column in \code{data} to
#'   plot. Probability summaries named \code{"Mean"}, \code{"Upper"}, or
#'   \code{"Lower"} are displayed on a fixed \eqn{[0, 1]} colour scale.
#' @param title Character string used as the plot title.
#' @param data Data frame containing columns \code{x1}, \code{x2}, and the
#'   summary variable named by \code{var}.
#' @param X Optional numeric matrix of observed two-dimensional input locations
#'   to overlay.
#' @param y Optional binary response vector corresponding to \code{X}. Points
#'   with \code{y == 1} and \code{y != 1} are drawn with different symbols.
#' @param dims Optional integer vector of length two indicating the original
#'   input dimensions being visualized. Used only for axis labels.
#' @param ... Additional arguments reserved for future use.
#'
#' @return A \code{ggplot} object.
#'
#' @keywords internal

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

#' Create a 2D ggplot surface for DKP class predictions or probabilities
#'
#' Internal helper used by \code{plot.DKP()} to create \pkg{ggplot2}-based
#' two-dimensional plots for predicted classes or maximum posterior class
#' probabilities. Observed class labels are overlaid using symbols determined by
#' the observed multinomial counts.
#'
#' @param var Character string giving the name of the column in \code{data} to
#'   plot. When \code{classification = TRUE}, the predicted class surface stored
#'   in \code{data$class} is used.
#' @param title Character string used as the plot title.
#' @param data Data frame containing columns \code{x1}, \code{x2}, and either a
#'   categorical \code{class} column or a continuous summary variable named by
#'   \code{var}.
#' @param X Numeric matrix of observed two-dimensional input locations to
#'   overlay.
#' @param Y Numeric matrix of multinomial counts or one-hot class indicators.
#'   The observed class is determined by \code{max.col(Y)}.
#' @param classification Logical. If \code{TRUE}, the plotted surface is treated
#'   as categorical class output; otherwise it is treated as a continuous
#'   probability or uncertainty surface.
#' @param dims Optional integer vector of length two indicating the original
#'   input dimensions being visualized. Used only for axis labels.
#' @param ... Additional arguments reserved for future use.
#'
#' @return A \code{ggplot} object.
#'
#' @keywords internal

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

#' Construct a fixed-slice plotting grid
#'
#' Internal helper used by plotting methods to expand a one- or two-dimensional
#' display grid into the model's full input dimension. Dimensions not displayed
#' in the plot are fixed at their training-data medians.
#'
#' @param X Training input matrix.
#' @param Xbounds Matrix of input bounds.
#' @param dims Integer vector of displayed dimensions.
#' @param grid Matrix or data frame containing grid values for displayed dimensions.
#' @return A numeric matrix with one row per grid point and one column per model input dimension.
#'
#' @keywords internal
.make_plot_grid <- function(X, Xbounds, dims, grid) {
  fixed <- apply(X, 2, stats::median)
  Xnew_full <- matrix(rep(fixed, each = nrow(grid)), nrow = nrow(grid))
  Xnew_full[, dims] <- as.matrix(grid)
  Xnew_full
}

#' Summarize posterior means and variances
#'
#' Internal helper used by summary methods to compute compact descriptive
#' statistics for posterior mean and variance vectors.
#'
#' @param mean_vals Numeric vector of posterior mean values.
#' @param var_vals Numeric vector of posterior variance values.
#'
#' @return A numeric matrix with rows for posterior means and posterior
#'   variances, and columns giving the mean, median, standard deviation, minimum,
#'   and maximum. Values are rounded to four decimal places.
#'
#' @keywords internal

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

#' Compute BKP posterior parameters
#'
#' Internal helper shared by \code{fit_BKP()}, \code{predict.BKP()}, and
#' \code{simulate.BKP()}. It computes kernel weights, prior Beta parameters,
#' posterior Beta parameters, and ESS diagnostics for a set of query locations.
#'
#' @param Xquery_norm Numeric matrix of normalized query locations.
#' @param Xtrain_norm Numeric matrix of normalized training locations.
#' @param y Numeric vector of observed success counts.
#' @param m Numeric vector of observed trial counts.
#' @param theta Positive kernel length-scale parameter or parameter vector.
#' @param kernel Character string specifying the kernel function.
#' @param isotropic Logical; whether to use a shared length-scale across input
#'   dimensions.
#' @param prior Character string specifying the prior type.
#' @param r0 Positive prior precision parameter.
#' @param p0 Prior mean parameter for fixed BKP priors.
#' @param ess Character string specifying the ESS calibration method.
#'
#' @return A list containing the kernel matrix, prior Beta parameters, posterior
#'   Beta parameters, and ESS diagnostics.
#'
#' @keywords internal

.bkp_compute_posterior <- function(Xquery_norm, Xtrain_norm, y, m, theta,
                                   kernel, isotropic, prior, r0, p0,
                                   ess = "none") {
  K <- kernel_matrix(Xquery_norm, Xtrain_norm, theta = theta,
                     kernel = kernel, isotropic = isotropic)
  prior_par <- get_prior(prior = prior, model = "BKP",
                         r0 = r0, p0 = p0, y = y, m = m, K = K)
  alpha0 <- prior_par$alpha0
  beta0 <- prior_par$beta0

  data_success <- as.vector(K %*% y)
  data_failure <- as.vector(K %*% (m - y))

  if (identical(ess, "shepard")) {
    ess_info <- .bkp_ess_calibration(
      Xquery_norm = Xquery_norm, Xtrain_norm = Xtrain_norm, m = m, K = K
    )
    data_success <- ess_info$scale * data_success
    data_failure <- ess_info$scale * data_failure
  } else {
    ess_info <- .bkp_ess_none_info(K, m)
  }

  list(
    K = K, alpha0 = alpha0, beta0 = beta0,
    alpha_n = alpha0 + data_success,
    beta_n = beta0 + data_failure,
    ess_info = ess_info
  )
}

#' Compute DKP posterior parameters
#'
#' Internal helper shared by \code{fit_DKP()}, \code{predict.DKP()}, and
#' \code{simulate.DKP()}. It computes kernel weights, prior Dirichlet
#' parameters, posterior Dirichlet parameters, and ESS diagnostics for a set of
#' query locations.
#'
#' @param Xquery_norm Numeric matrix of normalized query locations.
#' @param Xtrain_norm Numeric matrix of normalized training locations.
#' @param Y Numeric matrix of observed multinomial counts.
#' @param theta Positive kernel length-scale parameter or parameter vector.
#' @param kernel Character string specifying the kernel function.
#' @param isotropic Logical; whether to use a shared length-scale across input
#'   dimensions.
#' @param prior Character string specifying the prior type.
#' @param r0 Positive prior precision parameter.
#' @param p0 Prior mean vector for fixed DKP priors.
#' @param ess Character string specifying the ESS calibration method.
#'
#' @return A list containing the kernel matrix, prior Dirichlet parameters,
#'   posterior Dirichlet parameters, and ESS diagnostics.
#'
#' @keywords internal

.dkp_compute_posterior <- function(Xquery_norm, Xtrain_norm, Y, theta,
                                   kernel, isotropic, prior, r0, p0,
                                   ess = "none") {
  K <- kernel_matrix(Xquery_norm, Xtrain_norm, theta = theta,
                     kernel = kernel, isotropic = isotropic)
  alpha0 <- get_prior(prior = prior, model = "DKP", r0 = r0,
                      p0 = p0, Y = Y, K = K)

  data_counts <- as.matrix(K %*% Y)
  m <- rowSums(Y)
  if (identical(ess, "shepard")) {
    ess_info <- .bkp_ess_calibration(Xquery_norm, Xtrain_norm, m, K)
    data_counts <- sweep(data_counts, 1L, ess_info$scale, "*")
  } else {
    ess_info <- .bkp_ess_none_info(K, m)
  }

  list(K = K, alpha0 = alpha0, alpha_n = alpha0 + data_counts, ess_info = ess_info)
}

#' Check uniqueness of normalized input locations
#'
#' Internal validation helper for Shepard ESS calibration. Strict Shepard
#' interpolation is ambiguous when duplicated training locations carry different
#' trial sizes, so duplicated rows are rejected before calibration.
#'
#' @param Xnorm Numeric matrix of normalized input locations.
#'
#' @return Invisibly returns \code{TRUE} if all rows are unique; otherwise throws
#'   an error.
#'
#' @keywords internal

.bkp_check_unique_locations <- function(Xnorm) {
  if (anyDuplicated(as.data.frame(Xnorm)) > 0L) {
    stop(
      paste0(
        "ESS calibration with ess = 'shepard' requires unique input locations; ",
        "duplicated rows in 'X' are not supported because strict Shepard ",
        "interpolation is not well-defined for duplicates."
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Interpolate trial sizes using Shepard weights
#'
#' Internal helper used by optional Shepard ESS calibration. It interpolates
#' trial sizes from normalized training locations to normalized query locations
#' using inverse-distance weights, preserving exact-match values.
#'
#' @param Xquery_norm Numeric matrix of normalized query locations.
#' @param Xtrain_norm Numeric matrix of normalized training locations.
#' @param m Numeric vector of training trial sizes.
#' @param power Positive numeric Shepard inverse-distance power. The default is
#'   \code{2}.
#'
#' @return A numeric vector of interpolated trial sizes, one per query location.
#'
#' @keywords internal

.bkp_shepard_m <- function(Xquery_norm, Xtrain_norm, m, power = 2) {
  Xquery_norm <- as.matrix(Xquery_norm)
  Xtrain_norm <- as.matrix(Xtrain_norm)
  m <- as.numeric(m)

  if (ncol(Xquery_norm) != ncol(Xtrain_norm)) {
    stop("'Xquery_norm' and 'Xtrain_norm' must have the same number of columns.", call. = FALSE)
  }
  if (length(m) != nrow(Xtrain_norm)) {
    stop("'m' must have the same length as the number of training locations.", call. = FALSE)
  }
  if (!is.numeric(power) || length(power) != 1L || !is.finite(power) || power <= 0) {
    stop("'power' must be a positive finite scalar.", call. = FALSE)
  }

  as.vector(shepard_m_rcpp(Xquery_norm, Xtrain_norm, m, power))
}

#' Leave-one-out Shepard interpolation of trial sizes
#'
#' Internal helper used during hyperparameter tuning with \code{ess = "shepard"}.
#' For each training location, it interpolates the trial size from all other
#' training locations and excludes the observation's own contribution.
#'
#' @param Xnorm Numeric matrix of normalized training locations.
#' @param m Numeric vector of training trial sizes.
#' @param power Positive numeric Shepard inverse-distance power. The default is
#'   \code{2}.
#'
#' @return A numeric vector of leave-one-out interpolated trial sizes.
#'
#' @keywords internal
.bkp_shepard_m_loo <- function(Xnorm, m, power = 2) {
  Xnorm <- as.matrix(Xnorm)
  m <- as.numeric(m)

  .bkp_check_unique_locations(Xnorm)

  n <- nrow(Xnorm)
  if (n < 2L) {
    stop(
      "ESS calibration with ess = 'shepard' requires at least two input locations for leave-one-out calibration.",
      call. = FALSE
    )
  }
  if (length(m) != n) {
    stop("'m' must have the same length as the number of input locations.", call. = FALSE)
  }
  if (!is.numeric(power) || length(power) != 1L || !is.finite(power) || power <= 0) {
    stop("'power' must be a positive finite scalar.", call. = FALSE)
  }

  as.vector(shepard_m_loo_rcpp(Xnorm, m, power))
}

#' Compute Shepard ESS calibration diagnostics
#'
#' Internal helper that rescales kernel-weighted data contributions so their
#' effective trial size matches \eqn{\rho(\mathbf{x}) m_S(\mathbf{x})}, where
#' \eqn{m_S(\mathbf{x})} is a Shepard interpolation of trial sizes and
#' \eqn{\rho(\mathbf{x})} is the maximum kernel similarity to the training set.
#' The scaling preserves kernel-weighted empirical proportions and changes only
#' the data contribution.
#'
#' @param Xquery_norm Numeric matrix of normalized query locations.
#' @param Xtrain_norm Numeric matrix of normalized training locations.
#' @param m Numeric vector of training trial sizes.
#' @param K Numeric kernel matrix between query and training locations.
#'
#' @return A list containing \code{scale}, \code{m_kernel},
#'   \code{m_shepard}, \code{m_target}, and \code{rho}.
#'
#' @keywords internal

.bkp_ess_calibration <- function(Xquery_norm, Xtrain_norm, m, K) {
  .bkp_check_unique_locations(Xtrain_norm)

  m_kernel <- as.vector(K %*% as.numeric(m))
  m_shepard <- .bkp_shepard_m(Xquery_norm, Xtrain_norm, m, power = 2)
  rho <- apply(K, 1L, max)
  m_target <- rho * m_shepard

  scale <- numeric(length(m_kernel))
  positive_kernel_mass <- m_kernel > 0
  scale[positive_kernel_mass] <- m_target[positive_kernel_mass] / m_kernel[positive_kernel_mass]

  list(
    scale = scale,
    m_kernel = m_kernel,
    m_shepard = m_shepard,
    m_target = m_target,
    rho = rho
  )
}

#' Construct ESS diagnostics for the default uncalibrated update
#'
#' Internal helper used when \code{ess = "none"}. It returns the same diagnostic
#' structure as Shepard calibration, with unit scale factors and no Shepard
#' target values.
#'
#' @param K Numeric kernel matrix between query and training locations.
#' @param m Numeric vector of training trial sizes.
#'
#' @return A list containing unit \code{scale}, kernel-weighted trial sizes, and
#'   placeholder Shepard calibration fields.
#'
#' @keywords internal
.bkp_ess_none_info <- function(K, m) {
  list(
    scale = rep(1, nrow(K)),
    m_kernel = as.vector(K %*% as.numeric(m)),
    m_shepard = NULL,
    m_target = NULL,
    rho = apply(K, 1L, max)
  )
}


#' Compute TwinBKP posterior parameters using the C++ backend
#'
#' Internal helper used by \code{fit_TwinBKP()} and \code{predict.TwinBKP()} to
#' compute TwinBKP posterior Beta shape parameters using the uncalibrated
#' row-wise global-local pseudo-count aggregation in \code{twin_bkp_posterior_rcpp()}.
#'
#' @param Xquery_norm Numeric matrix of normalized query locations.
#' @param Xtrain_norm Numeric matrix of normalized training locations.
#' @param y Numeric vector of observed success counts.
#' @param m Numeric vector of observed binomial trial counts.
#' @param g_indices Integer vector of 1-based global subset indices.
#' @param local_indices Integer matrix of 1-based local-neighbour indices, with
#'   one row per query location.
#' @param theta_g Numeric scalar or vector giving the global kernel lengthscale.
#' @param theta_l Positive scalar giving the local kernel range.
#' @param global_kernel Character string specifying the global kernel.
#' @param local_kernel Character string specifying the local kernel.
#' @param isotropic Logical. If \code{TRUE}, the global kernel uses one shared
#'   lengthscale; otherwise it uses dimension-specific lengthscales.
#' @param prior Character string specifying the prior type.
#' @param r0 Positive scalar prior precision.
#' @param p0 Prior mean used by the fixed prior and passed to the C++ backend.
#' @param store_kernel Logical. If \code{TRUE}, dense diagnostic kernel matrices
#'   are returned; otherwise they are omitted.
#'
#' @return A list returned by \code{twin_bkp_posterior_rcpp()}, containing
#'   \code{alpha0}, \code{beta0}, \code{alpha_n}, \code{beta_n},
#'   \code{K}, \code{K_global}, and
#'   \code{K_local}.
#'
#' @keywords internal
.twin_bkp_compute_posterior <- function(
    Xquery_norm, Xtrain_norm, y, m, g_indices, local_indices,
    theta_g, theta_l, global_kernel, local_kernel, isotropic,
    prior, r0, p0, store_kernel = FALSE
) {
  twin_bkp_posterior_rcpp(
    Xquery_norm = Xquery_norm,
    Xtrain_norm = Xtrain_norm,
    y = as.numeric(y),
    m = as.numeric(m),
    g_indices = as.integer(g_indices),
    local_indices = local_indices,
    theta_g = as.numeric(theta_g),
    theta_l = as.numeric(theta_l),
    global_kernel = global_kernel,
    local_kernel = local_kernel,
    isotropic = isTRUE(isotropic),
    prior = prior,
    r0 = r0,
    p0 = p0,
    store_kernel = store_kernel
  )
}
