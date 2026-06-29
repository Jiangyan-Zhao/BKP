#' @name fit_TwinDKP
#'
#' @title Fit a Twin Dirichlet Kernel Process Model
#'
#' @description Fit a Twin Dirichlet Kernel Process (TwinDKP) model for
#'   categorical or multinomial count response data. TwinDKP is a scalable
#'   global-local approximation to the full \code{\link{fit_DKP}} model. It uses
#'   a twinning-selected global subset to capture the broad input-response
#'   structure and local nearest-neighbour updates to refine posterior inference
#'   near each target location.
#'
#' @inheritParams fit_DKP
#'
#' @param global_kernel Kernel function used for the global subset contribution.
#'   Options are \code{"gaussian"} (default), \code{"matern52"},
#'   \code{"matern32"}, and \code{"wendland"}.
#' @param local_kernel Kernel function used for the local-neighbour
#'   contribution. Currently only \code{"wendland"} is supported, corresponding
#'   to the compactly supported local kernel used by the TwinDKP approximation.
#' @param n_multi_start Number of initial points used in multi-start
#'   optimization of the global kernel lengthscale parameter(s). If
#'   \code{NULL}, the default from \code{\link{fit_DKP}} is used on the selected
#'   global subset.
#' @param n_threads Number of OpenMP threads used for global-subset
#'   hyperparameter optimization when \code{theta_g = NULL}. Default is
#'   \code{1}.
#' @param theta_g Optional positive global kernel lengthscale parameter. When
#'   \code{isotropic = TRUE}, this must be a scalar. When
#'   \code{isotropic = FALSE}, this can be either a scalar, which is broadcast
#'   to all dimensions, or a numeric vector of length \code{d}. If \code{NULL},
#'   the global lengthscale parameter is selected by fitting a DKP model on the
#'   selected global subset.
#' @param theta_l Optional positive scalar specifying the local Wendland kernel
#'   range. If \code{NULL}, it is set to the empirical covering radius of the
#'   global subset on the normalized input scale.
#' @param g Target global subset size. If \code{NULL}, the default is
#'   \eqn{\min\{n - 1, 50d, \max(\lfloor \sqrt{n} \rfloor, 10d)\}}.
#' @param l Number of local non-global neighbours used at each training or
#'   prediction location. If \code{NULL}, the default is
#'   \eqn{\min\{n - |G|, \max(25, 3d)\}} after the global subset has been
#'   selected.
#' @param twins Number of Twinning runs used to identify the global subset.
#'   Larger values may improve the selected global subset at additional
#'   computational cost. Default is \code{5}.
#' @param store_kernel Logical. If \code{TRUE}, store dense diagnostic kernel
#'   matrices \code{K}, \code{K_global}, and \code{K_local}. This option is
#'   intended for testing and diagnostics only. The default \code{FALSE} avoids
#'   dense \eqn{n \times n} kernel storage and preserves the scalable memory
#'   behavior of TwinDKP.
#'
#' @return A list of class \code{"TwinDKP"} containing the fitted TwinDKP model,
#'   with the following components:
#' \describe{
#'   \item{\code{theta_opt}}{Alias for \code{theta_g}, retained for consistency
#'     with \code{\link{fit_DKP}}.}
#'   \item{\code{theta_g}}{Optimized or user-specified global kernel lengthscale
#'     parameter(s).}
#'   \item{\code{theta_l}}{Local kernel range parameter used for the
#'     nearest-neighbour local component.}
#'   \item{\code{kernel}}{Alias for \code{global_kernel}, retained for
#'     consistency with \code{\link{fit_DKP}}.}
#'   \item{\code{global_kernel}}{Kernel function used for the global subset
#'     contribution.}
#'   \item{\code{local_kernel}}{Kernel function used for the local-neighbour
#'     contribution.}
#'   \item{\code{isotropic}}{Logical flag indicating whether the global kernel
#'     uses one shared lengthscale or per-dimension lengthscales.}
#'   \item{\code{loss}}{Loss function used for global lengthscale tuning.}
#'   \item{\code{loss_min}}{Loss value at the selected or user-specified global
#'     lengthscale parameter(s).}
#'
#'   \item{\code{X}}{Original training input matrix.}
#'   \item{\code{Xnorm}}{Training input matrix normalized to \eqn{[0,1]^d}.}
#'   \item{\code{Xbounds}}{Normalization bounds for each input dimension.}
#'   \item{\code{Y}}{Observed multinomial count matrix.}
#'
#'   \item{\code{prior}}{Prior specification used in the TwinDKP posterior
#'     update.}
#'   \item{\code{r0}}{Prior precision parameter.}
#'   \item{\code{p0}}{Prior class-probability vector used when
#'     \code{prior = "fixed"}.}
#'   \item{\code{alpha0}}{Prior Dirichlet concentration parameters evaluated at
#'     the training inputs.}
#'
#'   \item{\code{alpha_n}}{Posterior Dirichlet concentration parameters
#'     evaluated at the training inputs.}
#'   \item{\code{prob}}{Posterior mean class-probability matrix, computed as
#'     normalized rows of \code{alpha_n}.}
#'
#'   \item{\code{global_indices}}{One-based indices of the selected global
#'     subset.}
#'   \item{\code{control}}{A list of fitting controls and realized Twinning
#'     settings, including \code{g_target}, \code{g}, \code{l}, \code{r},
#'     \code{twins}, \code{u1}, \code{leaf_size}, \code{n_multi_start},
#'     \code{n_threads}, and \code{store_kernel}.}
#'   \item{\code{diagnostics}}{A list of diagnostic objects. It contains
#'     \code{twin_info} from the C++ Twinning routine and, when
#'     \code{store_kernel = TRUE}, dense matrices \code{K}, \code{K_global},
#'     and \code{K_local}. When \code{store_kernel = FALSE}, these kernel
#'     matrices are \code{NULL}.}
#' }
#'
#' @details TwinDKP first normalizes the input matrix to \eqn{[0,1]^d}. The
#'   global subset is selected by \code{twin_select_global_rcpp()} using the
#'   augmented representation \code{cbind(Xnorm, Y / rowSums(Y))}, so the
#'   selected points represent both the normalized input distribution and the
#'   empirical class-probability surface.
#'
#'   Given the selected global subset \eqn{G}, TwinDKP uses the union of
#'   \eqn{G} and a location-specific set of \code{l} nearest non-global
#'   neighbours for posterior aggregation. The global contribution uses
#'   \code{global_kernel} and \code{theta_g}; the local contribution uses the
#'   compactly supported \code{local_kernel} and \code{theta_l}. Local
#'   neighbours are found with a kd-tree over the non-global training points via
#'   \pkg{nanoflann}.
#'
#'   If \code{theta_g = NULL}, the global lengthscale parameter is selected by
#'   leave-one-out cross-validation on the selected global subset, using the
#'   specified \code{loss}. If \code{theta_g} is supplied, global tuning is
#'   skipped and the supplied value is used.
#'
#'   By default, TwinDKP aggregates posterior pseudo-counts row-wise and does
#'   not store dense \eqn{n \times n} kernel matrices. Fitting posterior
#'   aggregation costs \eqn{O(n(g + l))}. Prediction at \eqn{t} new input points
#'   costs \eqn{O(t(\log n + g + l))} for fixed input dimension. When
#'   \code{store_kernel = TRUE}, dense diagnostic matrices are additionally
#'   stored and memory use increases to \eqn{O(n^2)}.
#'
#'   Effective-sample-size calibration is currently available for full BKP and
#'   DKP models. TwinDKP uses the uncalibrated global-local posterior update to
#'   preserve the intended scalable approximation.
#'
#' @seealso \code{\link{fit_DKP}} for the full DKP model,
#'   \code{\link{fit_TwinBKP}} for the binomial TwinBKP analogue,
#'   \code{\link{fit_BKP}} for binary or binomial responses, and
#'   \code{\link{predict.TwinDKP}}, \code{\link{plot.TwinDKP}},
#'   \code{\link{simulate.TwinDKP}}, and \code{\link{summary.TwinDKP}} for
#'   downstream methods.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}. arXiv. <doi:10.48550/arXiv.2508.10447>.
#'
#' Vakayil A, Joseph VR (2022). Data Twinning. \emph{Statistical Analysis and
#'   Data Mining: The ASA Data Science Journal}, 15(5), 598--610.
#'   <doi:10.1002/sam.11574>.
#'
#' Vakayil A, Joseph VR (2024). A Global-Local Approximation Framework for
#'   Large-Scale Gaussian Process Modeling. \emph{Technometrics}, 66(2),
#'   295--305. <doi:10.1080/00401706.2023.2296451>.
#'
#' Blanco JL, PK Rai (2014). nanoflann: a C++ header-only fork of FLANN,
#'   a library for nearest neighbor (NN) with kd-trees.
#'   \url{https://github.com/jlblancoc/nanoflann}
#'
#' @examples
#' #-------------------------- 1D Example ---------------------------
#' set.seed(123)
#'
#' # Define true class probability function (3-class)
#' true_pi_fun <- function(X) {
#'   p1 <- 1/(1+exp(-3*X))
#'   p2 <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
#'   return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
#' }
#'
#' n <- 1000
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(150, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit TwinDKP model
#' model1 <- fit_TwinDKP(X, Y, Xbounds = Xbounds)
#' print(model1)
#'
#'
#' #-------------------------- 2D Example ---------------------------
#' # Define latent function and transform to 3-class probabilities
#' true_pi_fun <- function(X) {
#'   if (is.null(nrow(X))) X <- matrix(X, nrow = 1)
#'   m <- 8.6928; s <- 2.4269
#'   x1 <- 4 * X[,1] - 2
#'   x2 <- 4 * X[,2] - 2
#'   a <- 1 + (x1 + x2 + 1)^2 *
#'     (19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2)
#'   b <- 30 + (2*x1 - 3*x2)^2 *
#'     (18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2)
#'   f <- (log(a*b)- m)/s
#'   p1 <- pnorm(f) # Transform to probability
#'   p2 <- sin(pi * X[,1]) * sin(pi * X[,2])
#'   return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
#' }
#'
#' n <- 1000
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(150, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit TwinDKP model
#' model2 <- fit_TwinDKP(X, Y, Xbounds = Xbounds)
#' print(model2)
#'
#' @export
fit_TwinDKP <- function(
    X, Y, Xbounds = NULL,
    prior = c("noninformative", "fixed", "adaptive"),
    r0 = 2, p0 = NULL,
    global_kernel = c("gaussian", "matern52", "matern32", "wendland"),
    local_kernel = c("wendland"),
    loss = c("brier", "log_loss"),
    n_multi_start = NULL,
    isotropic = TRUE,
    n_threads = 1,
    theta_g = NULL,
    theta_l = NULL,
    g = NULL,
    l = NULL,
    twins = 5,
    store_kernel = FALSE
) {
  # ---- Argument checking ----
  if (missing(X) || missing(Y)) {
    stop("Arguments 'X' and 'Y' must be provided.")
  }
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("'X' must be a numeric matrix or data frame.")
  }
  if (!is.numeric(as.matrix(X))) {
    stop("'X' must contain numeric values only.")
  }
  if (!is.matrix(Y) && !is.data.frame(Y)) {
    stop("'Y' must be a numeric matrix or data frame.")
  }
  if (!is.numeric(as.matrix(Y))) {
    stop("'Y' must contain numeric values only.")
  }

  X <- as.matrix(X)
  Y <- as.matrix(Y)

  d <- ncol(X)
  n <- nrow(X)
  q <- ncol(Y)

  if (n < 3L) {
    stop("TwinDKP requires at least three observations.")
  }
  if (nrow(Y) != n) {
    stop("Number of rows in 'Y' must match number of rows in 'X'.")
  }
  if (q < 2L) {
    stop("'Y' must have at least two columns (multinomial outcomes).")
  }
  if (any(Y < 0)) {
    stop("'Y' must be nonnegative counts or frequencies.")
  }
  if (any(rowSums(Y) <= 0)) {
    stop("Each row of 'Y' must have a positive row sum.")
  }
  if (anyNA(X) || anyNA(Y)) {
    stop("Missing values are not allowed in 'X' or 'Y'.")
  }
  if (any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("'X' and 'Y' must contain only finite values.")
  }
  # ---- Prior, kernels, and loss ----
  prior <- match.arg(prior)
  global_kernel <- match.arg(global_kernel)
  local_kernel <- match.arg(local_kernel)
  loss <- match.arg(loss)

  if (is.null(p0)) {
    p0 <- colMeans(sweep(Y, 1, rowSums(Y), "/"))
  }

  if (!is.numeric(r0) || length(r0) != 1 || r0 <= 0 || is.na(r0) || !is.finite(r0)) {
    stop("'r0' must be a positive scalar.")
  }
  if (prior == "fixed" &&
      (is.null(p0) || !is.numeric(p0) || length(p0) != q || anyNA(p0) ||
       any(!is.finite(p0)) || any(p0 < 0) || abs(sum(p0) - 1) > 1e-10)) {
    stop("For fixed prior in DKP, 'p0' must be a nonnegative finite numeric vector of length equal to the number of classes and sum to 1.")
  }
  if (!is.logical(isotropic) || length(isotropic) != 1) {
    stop("'isotropic' must be a single logical value.")
  }

  if (!is.null(n_multi_start)) {
    if (!is.numeric(n_multi_start) || length(n_multi_start) != 1 ||
        is.na(n_multi_start) || !is.finite(n_multi_start) ||
        n_multi_start <= 0 || n_multi_start != floor(n_multi_start)) {
      stop("'n_multi_start' must be a positive integer.")
    }
    n_multi_start <- as.integer(n_multi_start)
  }

  if (!is.numeric(n_threads) || length(n_threads) != 1 ||
      is.na(n_threads) || !is.finite(n_threads) || n_threads <= 0 ||
      n_threads != floor(n_threads)) {
    stop("'n_threads' must be a positive integer.")
  }
  n_threads <- as.integer(n_threads)


  if (!is.null(theta_g)) {
    if (!is.numeric(theta_g)) {
      stop("'theta_g' must be numeric.")
    }
    if (isTRUE(isotropic) && length(theta_g) != 1) {
      stop("When isotropic=TRUE, 'theta_g' must be a scalar.")
    }
    if (!isTRUE(isotropic) && !(length(theta_g) == 1 || length(theta_g) == d)) {
      stop(paste0("When isotropic=FALSE, 'theta_g' must be either a scalar or a vector of length ", d, "."))
    }
    if (!isTRUE(isotropic) && length(theta_g) == 1) {
      theta_g <- rep(theta_g, d)
    }
    if (anyNA(theta_g) || any(!is.finite(theta_g)) || any(theta_g <= 0)) {
      stop("'theta_g' must be strictly positive.")
    }
  }
  if (!is.null(theta_l) &&
      (!is.numeric(theta_l) || length(theta_l) != 1 || is.na(theta_l) ||
       !is.finite(theta_l) || theta_l <= 0)) {
    stop("'theta_l' must be a positive scalar.")
  }
  if (!is.numeric(twins) || length(twins) != 1 || is.na(twins) ||
      !is.finite(twins) || twins <= 0 || twins != floor(twins)) {
    stop("'twins' must be a positive integer.")
  }
  if (!is.logical(store_kernel) || length(store_kernel) != 1) {
    stop("'store_kernel' must be a single logical value.")
  }

  # ---- Normalize input X to [0,1]^d ----
  if (is.null(Xbounds)) {
    xmin <- min(X)
    xmax <- max(X)
    if (xmin < 0 || xmax > 1) {
      warning(sprintf(
        "Input X does not appear to be normalized to [0,1]. Current range: [%.3f, %.3f].\nPlease normalize X or specify Xbounds explicitly; otherwise the model may produce incorrect results.",
        xmin, xmax
      ))
    }
    Xbounds <- cbind(rep(0, d), rep(1, d))
  } else {
    if (!is.matrix(Xbounds)) {
      stop("'Xbounds' must be a numeric matrix.")
    }
    if (!is.numeric(Xbounds)) {
      stop("'Xbounds' must contain numeric values.")
    }
    if (!all(dim(Xbounds) == c(d, 2))) {
      stop(paste0("'Xbounds' must be a matrix with dimensions d x 2, where d = ", d, "."))
    }
    if (anyNA(Xbounds) || any(!is.finite(Xbounds))) {
      stop("'Xbounds' must contain only finite values.")
    }
    if (any(Xbounds[, 2] <= Xbounds[, 1])) {
      stop("Each row of 'Xbounds' must satisfy lower < upper.")
    }
  }

  Xnorm <- sweep(X, 2, Xbounds[, 1], "-")
  Xnorm <- sweep(Xnorm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

  # ---- Twinning global subset selection ----
  if (is.null(g)) {
    g <- as.integer(max(2L, min(n - 1L, 50L * d, max(floor(sqrt(n)), 10L * d))))
  } else {
    if (!is.numeric(g) || length(g) != 1 || g != floor(g) ||
        is.na(g) || !is.finite(g) || g < 2 || g >= n) {
      stop("'g' must be an integer between 2 and n - 1.")
    }
    g <- as.integer(g)
  }

  r <- as.integer(max(2L, ceiling(n / g)))
  twins <- as.integer(twins)
  leaf_size <- 8L
  P_hat <- sweep(Y, 1, rowSums(Y), "/")
  twin_data <- cbind(Xnorm, P_hat)
  storage.mode(twin_data) <- "double"
  center <- colMeans(twin_data)
  first <- which.max(rowSums(sweep(twin_data, 2, center, "-")^2))
  u1 <- if (twins <= 1L) {
    as.integer(first)
  } else {
    as.integer(c(
      first,
      sample(
        setdiff(seq_len(n), first),
        twins - 1L,
        replace = (twins - 1L) > (n - 1L)
      )
    ))
  }

  twin_info <- twin_select_global_rcpp(
    twin_data = twin_data,
    Xnorm = Xnorm,
    r = r,
    runs = twins,
    u1 = u1,
    leaf_size = leaf_size
  )
  g_indices <- as.integer(twin_info$g_indices)
  g_actual <- length(g_indices)

  if (g_actual < 2L) {
    stop("The Twinning step selected fewer than two global points; try increasing 'g'.")
  }
  if (is.null(theta_l)) {
    theta_l <- as.numeric(twin_info$theta_l)
  }
  if (!is.finite(theta_l) || theta_l <= 0) {
    stop("The local range 'theta_l' must be positive.")
  }

  non_global_n <- n - g_actual
  if (is.null(l)) {
    l <- min(non_global_n, max(25L, 3L * d))
  } else {
    if (!is.numeric(l) || length(l) != 1 || l != floor(l) ||
        is.na(l) || !is.finite(l) || l < 0) {
      stop("'l' must be a nonnegative integer.")
    }
    l <- as.integer(l)
    if (l > non_global_n) {
      stop("'l' cannot exceed the number of non-global training points.")
    }
  }

  # ---- Global lengthscale tuning ----
  if (is.null(theta_g)) {
    global_fit <- fit_DKP(
      X[g_indices, , drop = FALSE],
      Y[g_indices, , drop = FALSE],
      Xbounds = Xbounds,
      prior = prior,
      r0 = r0,
      p0 = p0,
      kernel = global_kernel,
      loss = loss,
      n_multi_start = n_multi_start,
      theta = NULL,
      isotropic = isotropic,
      n_threads = n_threads,
      ess = "none"
    )
    theta_g <- as.numeric(global_fit$theta_opt)
    loss_min <- as.numeric(global_fit$loss_min)
  } else {
    theta_g <- as.numeric(theta_g)
    loss_min <- loss_fun(
      gamma = log10(theta_g),
      Xnorm = Xnorm[g_indices, , drop = FALSE],
      Y = Y[g_indices, , drop = FALSE],
      prior = prior,
      r0 = r0,
      p0 = p0,
      model = "DKP",
      loss = loss,
      kernel = global_kernel,
      isotropic = isotropic,
      ess = "none"
    )
  }

  local_indices <- twin_local_indices_rcpp(
    Xtrain_norm = Xnorm,
    Xquery_norm = Xnorm,
    g_indices = g_indices,
    l = l,
    leaf_size = leaf_size
  )

  # ---- Compute posterior parameters ----
  posterior <- twin_dkp_posterior_rcpp(
    Xquery_norm = Xnorm,
    Xtrain_norm = Xnorm,
    Y = Y,
    g_indices = as.integer(g_indices),
    local_indices = local_indices,
    theta_g = as.numeric(theta_g),
    theta_l = as.numeric(theta_l),
    global_kernel = global_kernel,
    local_kernel = local_kernel,
    isotropic = isTRUE(isotropic),
    prior = prior,
    r0 = r0,
    p0 = as.numeric(p0),
    store_kernel = store_kernel
  )

  # ---- Construct and return fitted model ----
  out <- list(
    # Kernel and tuning information
    theta_opt = theta_g,
    theta_g = theta_g,
    theta_l = theta_l,
    kernel = global_kernel,
    global_kernel = global_kernel,
    local_kernel = local_kernel,
    isotropic = isotropic,
    loss = loss,
    loss_min = loss_min,

    # Training data and normalization
    X = X,
    Xnorm = Xnorm,
    Xbounds = Xbounds,
    Y = Y,

    # Prior and posterior parameters
    prior = prior,
    r0 = r0,
    p0 = p0,
    alpha0 = posterior$alpha0,
    alpha_n = posterior$alpha_n,
    prob = posterior$prob,

    # Selected global subset
    global_indices = g_indices,

    # Fitting controls
    control = list(
      g_target = g,
      g = g_actual,
      l = l,
      r = r,
      twins = twins,
      u1 = u1,
      leaf_size = leaf_size,
      n_multi_start = n_multi_start,
      n_threads = n_threads,
      store_kernel = store_kernel
    ),

    # Diagnostic objects
    diagnostics = list(
      twin_info = twin_info,
      K = posterior$K,
      K_global = posterior$K_global,
      K_local = posterior$K_local
    )
  )
  class(out) <- "TwinDKP"
  out
}
