#' @name simulate
#'
#' @keywords DKP
#'
#' @examples
#' # -------------------------- DKP and TwinDKP ---------------------------
#' # Define true class probability function (3-class)
#' true_pi_fun <- function(X) {
#'   p1 <- 1/(1+exp(-3*X))
#'   p2 <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
#'   return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
#' }
#'
#' n <- 30
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(150, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit DKP model
#' model <- fit_DKP(X, Y, Xbounds = Xbounds, theta = 0.3)
#'
#' # Simulate 5 draws from posterior Dirichlet distributions at new point
#' Xnew <- matrix(seq(-2, 2, length.out = 5), ncol = 1)
#' simulate(model, Xnew = Xnew, nsim = 5)
#'
#' \dontrun{
#' # Larger TwinDKP example
#' n <- 200
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(150, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit TwinDKP model
#' model <- fit_TwinDKP(
#'      X, Y,
#'      Xbounds = Xbounds,
#'      theta_g = 0.3,
#'      g = 20,
#'      twins = 1,
#'      n_threads = 1
#'    )
#'
#' # Simulate 5 draws from posterior Dirichlet distributions at new point
#' simulate(model, Xnew = Xnew, nsim = 5)
#' }
#'
#' @export
#' @method simulate DKP

simulate.DKP <- function(object, nsim = 1, seed = NULL, Xnew = NULL, ...)
{
  # ---------------- Argument Checking ----------------
  if (!is.numeric(nsim) || length(nsim) != 1L ||
      is.na(nsim) || !is.finite(nsim) ||
      nsim <= 0 || nsim != as.integer(nsim)) {
    stop("`nsim` must be a positive integer.")
  }
  nsim <- as.integer(nsim)

  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L ||
       is.na(seed) || !is.finite(seed) ||
       seed != as.integer(seed))) {
    stop("`seed` must be a single integer or NULL.")
  }

  d <- ncol(object$X)

  if (!is.null(Xnew)) {
    if (is.null(dim(Xnew))) {
      if (d == 1L) {
        Xnew <- matrix(Xnew, ncol = 1L)
      } else {
        Xnew <- matrix(Xnew, nrow = 1L)
      }
    } else {
      Xnew <- as.matrix(Xnew)
    }

    if (!is.numeric(Xnew)) {
      stop("'Xnew' must be numeric.")
    }
    if (nrow(Xnew) < 1L || ncol(Xnew) < 1L) {
      stop("'Xnew' must have at least one row and one column.")
    }
    if (ncol(Xnew) != d) {
      stop("The number of columns in 'Xnew' must match the original input dimension.")
    }
    if (anyNA(Xnew) || any(!is.finite(Xnew))) {
      stop("'Xnew' must contain only finite values with no NA, NaN, or Inf.")
    }
  }

  # ---------------- Core Computation ----------------
  if (!is.null(seed)) set.seed(seed)

  if (!is.null(Xnew)) {
    prediction <- predict.DKP(object, Xnew = Xnew, type = "probability", ...)
    alpha_n <- prediction$alpha_n
    Xnew <- prediction$Xnew
    Y <- object$Y
    q <- ncol(Y)
    ess_info <- prediction$ess_info
  } else {
    # Use training data
    q       <- ncol(object$Y)
    alpha_n <- object$alpha_n
    Y       <- object$Y
    ess_info <- object$ess_info
  }


  # --- Simulate from Dirichlet posterior ---
  n_new <- ifelse(!is.null(Xnew), nrow(Xnew), nrow(object$X))
  samples <- array(0, dim = c(n_new, q, nsim))
  for (i in 1:n_new) {
    samples[i,,] <- t(rdirichlet(n = nsim, alpha = alpha_n[i, ]))  # [q × nsim]
  }

  dimnames(samples) <- list(
    paste0("x", 1:n_new),
    paste0("Class", 1:q),
    paste0("sim", 1:nsim)
  )

  # --- Optional: MAP prediction (only if data are single-label multinomial) ---
  class_pred <- NULL
  if (all(rowSums(Y) == 1)) {
    class_pred <- matrix(NA, nrow = n_new, ncol = nsim)
    for (i in 1:nsim) {
      class_pred[, i] <- max.col(samples[,,i])  # [n_new]
    }
  }

  # --- Posterior mean ---
  pi_mean <- alpha_n / rowSums(alpha_n)

  simulation <- list(
    samples = samples,    # [n_new × q × nsim]: posterior samples
    mean    = pi_mean,    # [n_new × q]: posterior mean
    class   = class_pred, # [n_new × nsim]: MAP class (if available)
    X       = object$X,   # [n × d]: training inputs
    Xnew    = Xnew,       # [n_new × d]: new inputs (if provided)
    ess     = if (is.null(object$ess)) "none" else object$ess,
    ess_info = ess_info
  )

  class(simulation) <- "simulate_DKP"
  return(simulation)
}
