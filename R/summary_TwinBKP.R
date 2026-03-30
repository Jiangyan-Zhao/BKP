#' @rdname summary
#' @keywords BKP
#'
#' @details
#' For \code{TwinBKP}, \code{summary()} reports global-stage information and
#' posterior summaries on the global support points selected by Twining.
#' This mirrors the role of \code{summary.BKP}, while respecting that TwinBKP
#' fits global hyperparameters first and performs full prediction behavior at
#' \code{predict.TwinBKP()} stage.
#'
#' @export
#' @method summary TwinBKP
summary.TwinBKP <- function(object, ...) {
  n_obs <- nrow(object$X)
  d <- ncol(object$X)
  gp <- .twinbkp_global_posterior(object)

  res <- list(
    n_obs = n_obs,
    input_dim = d,
    n_global = length(object$global_idx),
    g_nums = object$g_nums,
    kernel = object$kernel,
    isotropic = object$isotropic,
    theta_global = object$theta_global,
    loss = object$loss,
    loss_global = object$loss_global,
    prior = object$prior,
    r0 = object$r0,
    p0 = object$p0,
    global_idx = object$global_idx,
    post_mean_global = gp$mean,
    post_var_global = gp$variance
  )

  class(res) <- "summary_TwinBKP"
  res
}
