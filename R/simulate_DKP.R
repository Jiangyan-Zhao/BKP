#' @name simulate
#'
#' @keywords DKP
#'
#' @examples
#' ## -------------------- DKP --------------------
#' set.seed(123)
#'
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
#' model <- fit_DKP(X, Y, Xbounds = Xbounds)
#'
#' # Simulate 5 draws from posterior Dirichlet distributions at new point
#' Xnew <- matrix(seq(-2, 2, length.out = 5), ncol = 1)
#' simulate(model, Xnew = Xnew, nsim = 5)
#'
#' @export
#' @method simulate DKP

simulate.DKP <- function(object, nsim = 1, seed = NULL, Xnew = NULL, ...)
{
  # ---------------- Argument Checking ----------------
  if (!is.numeric(nsim) || length(nsim) != 1 || nsim <= 0 || nsim != as.integer(nsim)) {
    stop("`nsim` must be a positive integer.")
  }
  nsim <- as.integer(nsim)

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1 || seed != as.integer(seed))) {
    stop("`seed` must be a single integer or NULL.")
  }

  d <- ncol(object$Xnorm)
  if (!is.null(Xnew)) {
    if (is.null(nrow(Xnew))) {
      Xnew <- matrix(Xnew, nrow = 1)
    }
    Xnew <- as.matrix(Xnew)
    if (!is.numeric(Xnew)) {
      stop("'Xnew' must be numeric.")
    }
    if (ncol(Xnew) != d) {
      stop("The number of columns in 'Xnew' must match the original input dimension.")
    }
  }

  # ---------------- Core Computation ----------------
  if (!is.null(seed)) set.seed(seed)

  if (!is.null(Xnew)) {
    # complete posterior parameters at new inputs
    # Extract components
    Xnorm   <- object$Xnorm
    Y       <- object$Y
    theta   <- object$theta_opt
    kernel  <- object$kernel
    prior   <- object$prior
    r0      <- object$r0
    p0      <- object$p0
    Xbounds <- object$Xbounds
    q       <- ncol(Y)

    # --- Normalize new inputs ---
    Xnew_norm <- sweep(Xnew, 2, Xbounds[, 1], "-")
    Xnew_norm <- sweep(Xnew_norm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

    # --- Compute kernel matrix ---
    K <- kernel_matrix(Xnew_norm, Xnorm, theta = theta, kernel = kernel)

    # --- Get Dirichlet prior ---
    alpha0 <- get_prior(prior = prior, model = "DKP",
                        r0 = r0, p0 = p0, Y = Y, K = K)

    # --- Posterior Dirichlet parameters ---
    alpha_n <- as.matrix(alpha0) + as.matrix(K %*% Y)
    alpha_n <- pmax(alpha_n, 1e-10) # Avoid numerical issues
  }else{
    # Use training data
    q       <- ncol(object$Y)
    alpha_n <- pmax(object$alpha_n, 1e-10) # Avoid numerical issues
    Y       <- object$Y
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
    Xnew    = Xnew        # [n_new × d]: new inputs (if provided)
  )

  class(simulation) <- "simulate_DKP"
  return(simulation)
}
