#' @name summary
#'
#' @title Summarize Fitted BKP Package Models
#'
#' @description Provides a structured summary of fitted \code{BKP},
#'   \code{DKP}, \code{TwinBKP}, and \code{TwinDKP} model objects. The summary
#'   reports model dimensions, kernel settings, optimized hyperparameters,
#'   prior specification, loss information, approximation settings when
#'   applicable, and posterior summaries at the training input locations.
#'
#' @param object A fitted model object of class \code{"BKP"}, \code{"DKP"},
#'   \code{"TwinBKP"}, or \code{"TwinDKP"}, typically returned by
#'   \code{\link{fit_BKP}}, \code{\link{fit_DKP}},
#'   \code{\link{fit_TwinBKP}}, or \code{\link{fit_TwinDKP}}.
#' @param ... Additional arguments passed to the generic \code{summary} method
#'   (currently unused).
#'
#' @return A list containing model configuration, prior information, kernel
#'   hyperparameters, loss information, and posterior summaries at the training
#'   inputs.
#'
#'   For \code{BKP} and \code{TwinBKP}, posterior summaries are returned as
#'   vectors of posterior means and variances for the latent success
#'   probabilities.
#'
#'   For \code{DKP} and \code{TwinDKP}, posterior summaries are returned as
#'   matrices of marginal posterior means and variances for class probabilities.
#'
#'   For \code{TwinBKP} and \code{TwinDKP}, the list additionally includes
#'   global-local approximation fields such as \code{global_kernel},
#'   \code{local_kernel}, \code{theta_g}, \code{theta_l},
#'   \code{global_size}, \code{global_target}, \code{local_size}, and
#'   \code{twins}.
#'
#' @seealso \code{\link{fit_BKP}}, \code{\link{fit_DKP}},
#'   \code{\link{fit_TwinBKP}}, and \code{\link{fit_TwinDKP}} for model fitting;
#'   \code{\link{predict.BKP}}, \code{\link{predict.DKP}},
#'   \code{\link{predict.TwinBKP}}, and \code{\link{predict.TwinDKP}} for
#'   posterior prediction.
#'
#' @references Zhao J, Qing K, Xu J (2025). \emph{BKP: An R Package for Beta
#'   Kernel Process Modeling}. arXiv. <doi:10.48550/arXiv.2508.10447>.
#'
#' @keywords BKP
#'
#' @examples
#' # -------------------------- BKP and TwinBKP ---------------------------
#' set.seed(123)
#'
#' # Define true success probability function
#' true_pi_fun <- function(x) {
#'   (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' }
#'
#' n <- 30
#' Xbounds <- matrix(c(-2, 2), nrow = 1)
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit BKP model
#' model <- fit_BKP(X, y, m, Xbounds = Xbounds)
#' summary(model)
#'
#' \dontrun{
#' # Larger TwinBKP example
#' n <- 1000
#' X <- tgp::lhs(n = n, rect = Xbounds)
#' true_pi <- true_pi_fun(X)
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#'
#' # Fit TwinBKP model
#' model <- fit_TwinBKP(X, y, m, Xbounds = Xbounds)
#' summary(model)
#' }
#'
#' @export
#' @method summary BKP

summary.BKP <- function(object, ...) {
  # Extract information from the BKP object
  n_obs <- nrow(object$X)
  d     <- ncol(object$X)

  # Posterior summaries at training points
  post_mean <- as.vector(object$alpha_n / (object$alpha_n + object$beta_n))
  post_var  <- as.vector(post_mean * (1 - post_mean) / (object$alpha_n + object$beta_n + 1))

  res <- list(
    n_obs       = n_obs,
    input_dim   = d,
    kernel      = object$kernel,
    isotropic   = object$isotropic,
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
