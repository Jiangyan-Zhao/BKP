#' @title Loss Function for BKP and DKP Models
#'
#' @description Computes the loss used to fit the BKP (binary outcome) or DKP
#'   (multi-class outcome) model. Supports the Brier score (mean squared error)
#'   and negative log-loss (cross-entropy), under different prior
#'   specifications.
#'
#' @inheritParams fit.BKP
#' @inheritParams fit.DKP
#' @param gamma A numeric vector of log-transformed kernel hyperparameters.
#' @param Xnorm A numeric matrix of normalized inputs (each column scaled to
#'   \code{[0,1]}).
#' @param model_type A character string, either \code{"BKP"} (binary) or
#'   \code{"DKP"} (multi-class).
#'
#' @return A single numeric value representing the total loss (to be minimized):
#'   either the mean Brier score (squared error) or the mean negative log-loss
#'   (cross-entropy), depending on \code{loss}.
#'
#' @seealso \code{\link{fit.BKP}}, \code{\link{fit.DKP}},
#'   \code{\link{get_prior}}, \code{\link{kernel_matrix}}
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv.
#'   https://doi.org/10.48550/arXiv.2508.10447.
#'
#' @examples
#' ## Binary example (BKP)
#' set.seed(123)
#' n <- 10
#' Xnorm <- matrix(runif(2 * n), ncol = 2)
#' m <- rep(10, n)
#' y <- rbinom(n, size = m, prob = runif(n))
#' loss_fun(gamma = rep(0, 2), Xnorm = Xnorm, y = y, m = m, model_type = "BKP")
#'
#' ## Multi-class example (DKP)
#' set.seed(123)
#' n <- 10
#' q <- 3
#' Xnorm <- matrix(runif(2 * n), ncol = 2)
#' Y <- matrix(rmultinom(n, size = 10, prob = rep(1/q, q)), nrow = n, byrow = TRUE)
#' loss_fun(gamma = rep(0, 2), Xnorm = Xnorm, Y = Y, model_type = "DKP")
#'
#' @export

loss_fun <- function(
    gamma, Xnorm,
    y = NULL, m = NULL, Y = NULL,
    model_type = c("BKP", "DKP"),
    prior = c("noninformative", "fixed", "adaptive"), r0 = 2, p0 = 0.5,
    loss = c("brier", "log_loss"),
    kernel = c("gaussian", "matern52", "matern32"))
{
  # Match and validate the loss and kernel arguments
  model_type <- match.arg(model_type)
  prior <- match.arg(prior)
  loss <- match.arg(loss)
  kernel <- match.arg(kernel)

  # Convert gamma to kernel hyperparameters (theta = 10^gamma)
  theta <- 10^gamma

  # Compute kernel matrix using specified kernel and theta
  K <- kernel_matrix(Xnorm, theta = theta, kernel = kernel)
  diag(K) <- 0  # Leave-One-Out Cross-Validation (LOOCV)

  if (model_type == "BKP") {
    ## -------- Binary case (Beta-Binomial) --------
    # get the prior parameters: alpha0(x) and beta0(x)
    prior_par <- get_prior(
      prior = prior, r0 = r0, p0 = p0,
      y = y, m = m, K = K
    )
    alpha0 <- prior_par$alpha0
    beta0 <- prior_par$beta0

    # Compute posterior alpha and beta
    alpha_n <- alpha0 + as.vector(K %*% y)
    beta_n <- beta0 + as.vector(K %*% (m - y))

    # Posterior mean prediction of success probability
    pi_hat <- alpha_n / (alpha_n + beta_n)
    pi_hat <- pmin(pmax(pi_hat, 1e-10), 1 - 1e-10)   # avoid log(0)

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
    # get the prior parameters: alpha0(x) and beta0(x)
    alpha0 <- get_prior(prior = prior, model_type = "DKP",
                        r0 = r0, p0 = p0, Y = Y, K = K)

    # Compute posterior alpha
    alpha_n <- as.matrix(alpha0) + as.matrix(K %*% Y)

    # Posterior mean prediction of success probability
    pi_hat <- alpha_n / rowSums(alpha_n)
    pi_hat <- pmin(pmax(pi_hat, 1e-6), 1 - 1e-6)   # avoid log(0)

    if (loss == "brier") {
      # Brier score (Mean Squared Error)
      # Empirical success rate
      pi_tilde <- Y / rowSums(Y)
      # Brier score: mean squared error between predicted and observed
      brier <- mean((pi_hat - pi_tilde)^2)
      return(brier)
    } else {
      # log-loss (cross-entropy)
      log_loss <- -mean(Y * log(pi_hat))
      return(log_loss)
    }
  }
}
