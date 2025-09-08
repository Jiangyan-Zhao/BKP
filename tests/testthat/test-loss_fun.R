#' @title Loss Function for BKP and DKP Models
#'
#' @description Computes the loss for fitting BKP (binary) or DKP (multi-class)
#'   models. Supports Brier score (mean squared error) and log-loss
#'   (cross-entropy) under different prior specifications.
#'
#' @inheritParams get_prior
#' @inheritParams fit_BKP
#' @param gamma A numeric vector of log10-transformed kernel hyperparameters.
#' @param Xnorm A numeric matrix of normalized input features ([0,1]^d).
#' @param model A character string specifying the model type: "BKP"
#'   (binary) or "DKP" (multi-class).
#'
#' @return A single numeric value representing the total loss (to be minimized).
#'   The value corresponds to either the Brier score (squared error) or the
#'   log-loss (cross-entropy).
#'
#' @seealso \code{\link{fit_BKP}} for fitting BKP models, \code{\link{fit_DKP}}
#'   for fitting DKP models, \code{\link{get_prior}} for constructing prior
#'   parameters, \code{\link{kernel_matrix}} for computing kernel matrices.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv.
#'   https://doi.org/10.48550/arXiv.2508.10447.
#'
#' @examples
#' # -------------------------- BKP ---------------------------
#' set.seed(123)
#' d <- 2
#' n <- 50
#' Xnorm <- matrix(runif(n * d), n, d)
#' y <- rbinom(n, size = 10, prob = 0.5)
#' m <- rep(10, n)
#' theta <- 0.5
#'
#' # Compute loss for a noninformative prior
#' loss_noninformative <- loss_fun(gamma = log10(theta), Xnorm = Xnorm, y = y, m = m,
#'                                 model = "BKP", loss = "brier", prior = "noninformative")
#'
#' # Compute loss for a fixed prior
#' loss_fixed <- loss_fun(gamma = log10(theta), Xnorm = Xnorm, y = y, m = m,
#'                        model = "BKP", loss = "brier", prior = "fixed", r0 = 10, p0 = 0.6)
#'
#' # -------------------------- DKP ---------------------------
#' set.seed(123)
#' d <- 2
#' n <- 50
#' q <- 3
#' Xnorm <- matrix(runif(n * d), n, d)
#' Y <- matrix(c(rbinom(n, size = 10, prob = 0.3),
#'               rbinom(n, size = 10, prob = 0.4),
#'               rbinom(n, size = 10, prob = 0.3)), n, q)
#'
#' # Compute loss for an adaptive prior
#' loss_adaptive <- loss_fun(gamma = log10(theta), Xnorm = Xnorm, Y = Y,
#'                           model = "DKP", loss = "brier", prior = "adaptive", r0 = 10)
#'
#' @export
loss_fun <- function(gamma, Xnorm, y = NULL, m = NULL, Y = NULL,
                     model = c("BKP", "DKP"), loss = c("brier", "log_loss"),
                     prior = c("noninformative", "fixed", "adaptive"),
                     r0 = NULL, p0 = NULL, kernel = "gaussian") {

  # Match arguments to ensure validity
  model <- match.arg(model)
  loss <- match.arg(loss)
  prior <- match.arg(prior)

  # Check inputs
  if (!is.null(y) && !is.null(m)) {
    if (length(y) != length(m)) {
      stop("'y' and 'm' must have the same length.")
    }
  }

  if (length(gamma) != ncol(Xnorm)) {
    stop("length(gamma) must be equal to ncol(Xnorm).")
  }

  theta <- 10^gamma

  # Compute kernel matrix K
  K <- kernel_matrix(Xnorm, theta = theta, kernel = kernel)

  # Set diagonal of K to zero for leave-one-out cross-validation
  diag(K) <- 0

  if (model == "BKP") {
    ## -------- Binary case (Beta-Binomial) --------
    # get the prior parameters: alpha0(x) and beta0(x)
    prior_par <- get_prior(prior = prior, model = model, r0 = r0, p0 = p0, y = y, m = m, K = K)
    alpha0 <- prior_par$alpha0
    beta0 <- prior_par$beta0

    # Compute posterior alpha_n and beta_n
    alpha_n <- alpha0 + as.vector(K %*% y)
    beta_n <- beta0 + as.vector(K %*% (m - y))

    # Posterior mean prediction of success probability
    pi_hat <- alpha_n / (alpha_n + beta_n)
    pi_hat <- pmin(pmax(pi_hat, 1e-10), 1 - 1e-10) # avoid log(0)

    if (loss == "brier") {
      # Brier score (Mean Squared Error)
      # Empirical success rate
      pi_tilde <- y / m
      # Brier score: mean squared error between predicted and observed
      brier <- mean((pi_hat - pi_tilde)^2)
      return(brier)
    } else {
      # log-loss (cross-entropy)
      log_loss <- -mean(y * log(pi_hat) + (m - y) * log(1 - pi_hat))
      return(log_loss)
    }
  } else {
    ## -------- Multiclass case (Dirichlet-Multinomial) --------
    # get the prior parameters: alpha0(x)
    alpha0 <- get_prior(prior = prior, model = model, r0 = r0, p0 = p0, Y = Y, K = K)

    # Compute posterior alpha
    alpha_n <- as.matrix(alpha0) + as.matrix(K %*% Y)

    # Posterior mean prediction of success probability
    pi_hat <- alpha_n / rowSums(alpha_n)
    pi_hat <- pmin(pmax(pi_hat, 1e-10), 1 - 1e-10) # avoid log(0)

    if (loss == "brier") {
      # Brier score (Mean Squared Error)
      Y_normalized <- Y / rowSums(Y)
      brier <- mean((pi_hat - Y_normalized)^2)
      return(brier)
    } else {
      # log-loss (cross-entropy)
      log_loss <- -mean(rowSums(Y * log(pi_hat)))
      return(log_loss)
    }
  }
}
