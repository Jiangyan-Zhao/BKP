#' @name simulate
#'
#' @keywords DKP
#'
#' @examples
#' ## -------------------- DKP Simulation Example --------------------
#' set.seed(123)
#'
#' # Define true class probability function (3-class)
#' true_pi_fun <- function(X) {
#'   p <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
#'   return(matrix(c(p/2, p/2, 1 - p), nrow = length(p)))
#' }
#'
#' n <- 30
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit DKP model
#' model <- fit.DKP(X, Y, Xbounds = Xbounds)
#'
#' # Simulate 5 draws from posterior Dirichlet distributions at new point
#' Xnew <- matrix(seq(-2, 2, length.out = 100), ncol = 1)
#' simulate(model, Xnew = Xnew, n_sim = 5)
#'
#' @export
#' @method simulate DKP

simulate.DKP <- function(object, Xnew, n_sim = 1, seed = NULL, ...) {
  if (!inherits(object, "DKP")) {
    stop("The input must be of class 'DKP'. Please provide a model fitted with 'fit.DKP()'.")
  }
  DKPmodel <- object

  # Extract components
  Xnorm   <- DKPmodel$Xnorm
  Y       <- DKPmodel$Y
  theta   <- DKPmodel$theta_opt
  kernel  <- DKPmodel$kernel
  prior   <- DKPmodel$prior
  r0      <- DKPmodel$r0
  p0      <- DKPmodel$p0
  Xbounds <- DKPmodel$Xbounds
  d       <- ncol(Xnorm)
  q       <- ncol(Y)

  # --- Check and normalize Xnew ---
  if (is.null(nrow(Xnew))) Xnew <- matrix(Xnew, nrow = 1)
  Xnew <- as.matrix(Xnew)
  if (ncol(Xnew) != d) {
    stop("Xnew must have the same number of columns as the training input.")
  }
  n_new <- nrow(Xnew)
  Xnew_norm <- sweep(Xnew, 2, Xbounds[, 1], "-")
  Xnew_norm <- sweep(Xnew_norm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

  # --- Compute kernel matrix ---
  K <- kernel_matrix(Xnew_norm, Xnorm, theta = theta, kernel = kernel)

  # --- Get Dirichlet prior ---
  alpha0 <- get_prior_dkp(prior = prior, r0 = r0, p0 = p0, Y = Y, K = K)

  # --- Posterior Dirichlet parameters ---
  alpha_n <- as.matrix(alpha0) + as.matrix(K %*% Y)
  alpha_n <- pmax(alpha_n, 1e-10) # Avoid numerical issues

  # --- Simulate from Dirichlet posterior ---
  if (!is.null(seed)) set.seed(seed)
  sims <- array(0, dim = c(n_sim, q, n_new))
  for (i in 1:n_new) {
    shape_mat <- matrix(rgamma(n_sim * q, shape = rep(alpha_n[i, ], each = n_sim), rate = 1),
                        nrow = n_sim)
    sims[,,i] <- shape_mat / rowSums(shape_mat)
  }

  dimnames(sims) <- list(
    paste0("sim", 1:n_sim),
    paste0("Class", 1:q),
    paste0("x", 1:n_new)
  )

  # --- Posterior mean ---
  pi_mean <- alpha_n / rowSums(alpha_n)

  # --- Optional: MAP prediction (only if data are single-label multinomial) ---
  class_pred <- NULL
  if (all(rowSums(Y) == 1)) {
    class_pred <- matrix(NA, nrow = n_new, ncol = n_sim)
    for (i in 1:n_sim) {
      class_pred[, i] <- max.col(t(sims[i,,]))  # [n_new]
    }
  }

  return(list(
    sims = sims,        # [n_sim × q × n_new]: posterior samples
    mean = pi_mean,     # [n_new × q]: posterior mean
    class = class_pred  # [n_new × n_sim]: MAP class (if available)
  ))
}

