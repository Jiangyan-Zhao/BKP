#' @rdname fitted
#' @keywords TwinDKP
#' @export
#' @method fitted TwinDKP
fitted.TwinDKP <- function(object, ...) {
  # assume n x q matrix (q classes)
  alpha_n <- object$alpha_n
  q <- ncol(alpha_n)

  # Posterior mean
  fitted_value <- alpha_n / rowSums(alpha_n)
  colnames(fitted_value) <- paste0("class", seq_len(q))

  return(fitted_value)
}
