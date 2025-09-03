#' @name summary
#'
#' @title Summary of a Fitted BKP or DKP Model
#'
#' @description Provides a structured summary of a fitted Beta Kernel Process
#'   (BKP) or Dirichlet Kernel Process (DKP) model. This function reports the
#'   model configuration, prior specification, kernel settings, and key
#'   posterior quantities, giving users a concise overview of the fitting
#'   results.
#'
#' @param object An object of class \code{"BKP"} (from \code{\link{fit_BKP}}) or
#'   \code{"DKP"} (from \code{\link{fit_DKP}}).
#' @param ... Additional arguments passed to the generic \code{summary} method
#'   (currently not used).
#'
#' @return A list containing key summaries of the fitted model:
#' \describe{
#'   \item{\code{n_obs}}{Number of training observations.}
#'   \item{\code{input_dim}}{Input dimensionality (number of columns in X).}
#'   \item{\code{kernel}}{Kernel type used in the model.}
#'   \item{\code{theta_opt}}{Estimated kernel hyperparameters.}
#'   \item{\code{loss}}{Loss function type used in the model.}
#'   \item{\code{loss_min}}{Minimum value of the loss function achieved.}
#'   \item{\code{prior}}{Prior type used (e.g., "noninformative", "fixed", "adaptive").}
#'   \item{\code{r0}}{Prior precision parameter.}
#'   \item{\code{p0}}{Prior mean parameter.}
#'   \item{\code{post_mean}}{Posterior mean estimates.
#'     For BKP: posterior mean success probabilities at training points.
#'     For DKP: posterior mean class probabilities (\eqn{n_\text{obs} \times q}).}
#'   \item{\code{post_var}}{Posterior variance estimates.
#'     For BKP: variance of success probabilities.
#'     For DKP: variance for each class probability.}
#'   \item{\code{n_class}}{(Only for DKP) Number of classes in the response.}
#' }
#'
#' @seealso \code{\link{fit_BKP}}, \code{\link{fit_DKP}} for model fitting.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}.  arXiv.
#'   https://doi.org/10.48550/arXiv.2508.10447.
#'
#' @keywords BKP
#'
#' @examples
#' # ============================================================== #
#' # ========================= BKP Examples ======================= #
#' # ============================================================== #
#'
#' #-------------------------- 1D Example ---------------------------
#' set.seed(123)
#'
#' # Define true success probability function
#' true_pi_fun <- function(x) {
#'   (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' }
#'
#' n <- 30
#' Xbounds <- matrix(c(-2,2), nrow=1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit BKP model
#' model1 <- fit_BKP(X, y, m, Xbounds=Xbounds)
#' summary(model1)
#'
#'
#' #-------------------------- 2D Example ---------------------------
#' set.seed(123)
#'
#' # Define 2D latent function and probability transformation
#' true_pi_fun <- function(X) {
#'   if(is.null(nrow(X))) X <- matrix(X, nrow=1)
#'   m <- 8.6928
#'   s <- 2.4269
#'   x1 <- 4*X[,1]- 2
#'   x2 <- 4*X[,2]- 2
#'   a <- 1 + (x1 + x2 + 1)^2 *
#'     (19- 14*x1 + 3*x1^2- 14*x2 + 6*x1*x2 + 3*x2^2)
#'   b <- 30 + (2*x1- 3*x2)^2 *
#'     (18- 32*x1 + 12*x1^2 + 48*x2- 36*x1*x2 + 27*x2^2)
#'   f <- log(a*b)
#'   f <- (f- m)/s
#'   return(pnorm(f))  # Transform to probability
#' }
#'
#' n <- 100
#' Xbounds <- matrix(c(0, 0, 1, 1), nrow = 2)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit BKP model
#' model2 <- fit_BKP(X, y, m, Xbounds=Xbounds)
#' summary(model2)
#'
#' @export
#' @method summary BKP

summary.BKP <- function(object, ...) {
  # Extract information from the BKP object
  n_obs <- nrow(object$X)
  d     <- ncol(object$X)

  # Posterior summaries at training points
  post_mean <- object$alpha_n / (object$alpha_n + object$beta_n)
  post_var  <- post_mean * (1 - post_mean) / (object$alpha_n + object$beta_n + 1)

  res <- list(
    n_obs       = n_obs,
    input_dim   = d,
    kernel      = object$kernel,
    theta_opt   = object$theta_opt,
    loss        = object$loss,
    loss_min    = object$loss_min,
    prior       = object$prior,
    r0          = object$r0,
    p0          = object$p0,
    post_mean   = post_mean,
    post_var    = post_var
  )

  class(res) <- "summary_BKP"
  return(res)
}
