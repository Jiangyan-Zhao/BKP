#' Default optimization control for BKP and DKP
#'
#' Constructs and validates default optimization controls for kernel
#' hyperparameter tuning. Different defaults are used for isotropic and
#' anisotropic kernels.
#'
#' @param d Input dimension.
#' @param isotropic Logical; whether an isotropic kernel is used.
#' @param control A user-supplied list of optimization controls.
#'
#' @return A validated list with elements \code{n_grid}, \code{n_starts},
#'   \code{max_iter}, \code{g_lower}, and \code{g_upper}.
#'
#' @keywords internal
.default_opt_control <- function(d, isotropic, control = list()) {
  # ---- Argument checking ----
  if (!is.numeric(d) || length(d) != 1L || is.na(d) || !is.finite(d)) {
    stop("'d' must be a finite positive integer.", call. = FALSE)
  }

  if (d < 1L || d != as.integer(d)) {
    stop("'d' must be a finite positive integer.", call. = FALSE)
  }

  d <- as.integer(d)

  if (!is.logical(isotropic) || length(isotropic) != 1L || is.na(isotropic)) {
    stop("'isotropic' must be TRUE or FALSE.", call. = FALSE)
  }

  if (is.null(control)) {
    control <- list()
  }

  if (!is.list(control)) {
    stop("'control' must be a list.", call. = FALSE)
  }

  # ---- Defaults ----
  if (isTRUE(isotropic)) {
    defaults <- list(
      n_grid = 100L,
      n_starts = 1L,
      max_iter = 300L,
      g_lower = -3,
      g_upper = 3
    )
  } else {
    defaults <- list(
      n_grid = as.integer(min(1000L, 100L * d)),
      n_starts = as.integer(min(100L, 10L * d)),
      max_iter = as.integer(min(500L, ceiling(100 * log1p(d)))),
      g_lower = -3,
      g_upper = 3
    )
  }

  # ---- Check control names ----
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown) > 0L) {
    stop(
      "Unknown control parameter(s): ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  out <- utils::modifyList(defaults, control)

  # ---- Validate control values ----
  if (!is.numeric(out$n_grid) || length(out$n_grid) != 1L ||
      is.na(out$n_grid) || !is.finite(out$n_grid)) {
    stop("'control$n_grid' must be a finite positive integer.", call. = FALSE)
  }

  if (!is.numeric(out$n_starts) || length(out$n_starts) != 1L ||
      is.na(out$n_starts) || !is.finite(out$n_starts)) {
    stop("'control$n_starts' must be a finite positive integer.", call. = FALSE)
  }

  if (!is.numeric(out$max_iter) || length(out$max_iter) != 1L ||
      is.na(out$max_iter) || !is.finite(out$max_iter)) {
    stop("'control$max_iter' must be a finite positive integer.", call. = FALSE)
  }

  out$n_grid <- as.integer(out$n_grid)
  out$n_starts <- as.integer(out$n_starts)
  out$max_iter <- as.integer(out$max_iter)

  if (out$n_grid < 5L) {
    stop("'control$n_grid' must be at least 5.", call. = FALSE)
  }

  if (out$n_starts < 1L) {
    stop("'control$n_starts' must be at least 1.", call. = FALSE)
  }

  if (out$max_iter < 1L) {
    stop("'control$max_iter' must be positive.", call. = FALSE)
  }

  if (out$n_starts > out$n_grid) {
    out$n_starts <- out$n_grid
  }

  if (!is.numeric(out$g_lower) || length(out$g_lower) != 1L ||
      is.na(out$g_lower) || !is.finite(out$g_lower)) {
    stop("'control$g_lower' must be a finite numeric scalar.", call. = FALSE)
  }

  if (!is.numeric(out$g_upper) || length(out$g_upper) != 1L ||
      is.na(out$g_upper) || !is.finite(out$g_upper)) {
    stop("'control$g_upper' must be a finite numeric scalar.", call. = FALSE)
  }

  if (out$g_lower >= out$g_upper) {
    stop(
      "'control$g_lower' must be smaller than 'control$g_upper'.",
      call. = FALSE
    )
  }

  out
}
