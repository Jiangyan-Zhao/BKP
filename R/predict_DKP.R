#' @rdname predict
#'
#' @keywords DKP
#'
#' @examples
#' # ============================================================== #
#' # ========================= DKP Examples ======================= #
#' # ============================================================== #
#'
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
#' model1 <- fit_DKP(X, Y, Xbounds = Xbounds)
#'
#' # Prediction on training data
#' predict(model1)
#'
#' # Prediction on new data
#' Xnew = matrix(seq(-2, 2, length = 10), ncol=1) #new data points
#' predict(model1, Xnew)
#'
#'
#' #-------------------------- 2D Example ---------------------------
#' set.seed(123)
#'
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
#' n <- 100
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(150, n, replace = TRUE)
#'
#' # Generate multinomial responses
#' Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))
#'
#' # Fit DKP model
#' model2 <- fit_DKP(X, Y, Xbounds = Xbounds)
#'
#' # Prediction on training data
#' predict(model2)
#'
#' # Prediction on new data
#' x1 <- seq(Xbounds[1,1], Xbounds[1,2], length.out = 10)
#' x2 <- seq(Xbounds[2,1], Xbounds[2,2], length.out = 10)
#' Xnew <- expand.grid(x1 = x1, x2 = x2)
#' predict(model2, Xnew)
#'
#' @export
#' @method predict DKP

predict.DKP <- function(object, Xnew = NULL, CI_level = 0.95,
                        return_type = c("probability", "count"), n_trials = NULL, ...)
{
  #---- Argument checks ----
  X       <- object$X
  d       <- ncol(X)
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

  if (!is.numeric(CI_level) || length(CI_level) != 1 || CI_level <= 0 || CI_level >= 1) {
    stop("'CI_level' must be a single numeric value strictly between 0 and 1.")
  }
  return_type <- match.arg(return_type)
  if (return_type == "count") {
    if (is.null(n_trials)) {
      stop("When return_type = 'count', 'n_trials' must be provided.")
    }
    if (!is.numeric(n_trials) || length(n_trials) != 1 || n_trials <= 0 || n_trials != as.integer(n_trials)) {
      stop("'n_trials' must be a positive integer when return_type = 'count'.")
    }
    n_trials <- as.integer(n_trials)
  }

  if(!is.null(Xnew)){
    # Extract components
    Xnorm   <- object$Xnorm
    Y       <- object$Y
    theta   <- object$theta_opt
    kernel  <- object$kernel
    isotropic <- object$isotropic
    prior   <- object$prior
    r0      <- object$r0
    p0      <- object$p0
    Xbounds <- object$Xbounds

    # Normalize Xnew to [0,1]^d
    Xnew_norm <- sweep(Xnew, 2, Xbounds[, 1], "-")
    Xnew_norm <- sweep(Xnew_norm, 2, Xbounds[, 2] - Xbounds[, 1], "/")

    # Compute kernel matrix
    K <- kernel_matrix(Xnew_norm, Xnorm, theta = theta, kernel = kernel, isotropic = isotropic)

    # Get prior parameters
    alpha0 <- get_prior(prior = prior, model = "DKP",
                        r0 = r0, p0 = p0, Y = Y, K = K)

    # Call C++ function for posterior computation
    result <- predict_dkp_rcpp(K, as.matrix(alpha0), as.matrix(Y))
    alpha_n <- result$alpha_n
    row_sum <- result$row_sum
    pred_mean <- result$mean
    pred_var  <- result$variance
  }else{
    # Use training data
    alpha_n <- object$alpha_n
    Y       <- object$Y
    row_sum <- rowSums(alpha_n)
    
    eps <- 1e-10
    pred_mean <- alpha_n / pmax(row_sum, eps)
    pred_mean <- pmin(pmax(pred_mean, eps), 1 - eps)
    pred_var  <- pred_mean * (1 - pred_mean) / (row_sum + 1)
  }

  beta_n <- matrix(row_sum, nrow = nrow(alpha_n), ncol = ncol(alpha_n)) - alpha_n
  if (return_type == "probability") {
    # Credible intervals on class probability p_j
    pred_lower <- suppressWarnings(qbeta((1 - CI_level) / 2, alpha_n, beta_n))
    pred_upper <- suppressWarnings(qbeta((1 + CI_level) / 2, alpha_n, beta_n))
  } else {
    # Count prediction via marginal Beta-Binomial(n_trials, alpha_j, row_sum-alpha_j)
    pred_mean <- n_trials * alpha_n / pmax(row_sum, 1e-10)
    pred_var <- n_trials * alpha_n * beta_n * (matrix(row_sum, nrow(alpha_n), ncol(alpha_n)) + n_trials) /
      (pmax(matrix(row_sum^2, nrow(alpha_n), ncol(alpha_n)), 1e-10) *
         pmax(matrix(row_sum + 1, nrow(alpha_n), ncol(alpha_n)), 1e-10))

    pred_lower <- matrix(0, nrow = nrow(alpha_n), ncol = ncol(alpha_n))
    pred_upper <- matrix(0, nrow = nrow(alpha_n), ncol = ncol(alpha_n))
    for (i in seq_len(nrow(alpha_n))) {
      for (j in seq_len(ncol(alpha_n))) {
        pred_lower[i, j] <- betabinom_quantile((1 - CI_level) / 2, n_trials, alpha_n[i, j], beta_n[i, j])
        pred_upper[i, j] <- betabinom_quantile((1 + CI_level) / 2, n_trials, alpha_n[i, j], beta_n[i, j])
      }
    }
  }
  class_names <- paste0("class", seq_len(ncol(alpha_n)))
  colnames(pred_mean) <- class_names
  colnames(pred_var) <- class_names
  colnames(pred_lower) <- class_names
  colnames(pred_upper) <- class_names

  # Return structured output
  prediction <- list(
    X        = X,
    Xnew     = Xnew,
    alpha_n  = alpha_n,
    mean     = pred_mean,
    variance = pred_var,
    lower    = pred_lower,
    upper    = pred_upper,
    CI_level = CI_level,
    return_type = return_type
  )

  if (return_type == "count") {
    prediction$n_trials <- n_trials
  }

  # Posterior classification label (only for classification data)
  if (return_type == "probability" && all(rowSums(Y) == 1)) {
    prediction$class <- max.col(pred_mean)
  }

  class(prediction) <- "predict_DKP"
  return(prediction)
}