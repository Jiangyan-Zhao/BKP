#' Internal Twin index selector
#' @keywords internal
get_twin_indices <- function(Xnorm, g, v = 2L * g, runs = 10L, seed = 123L) {
  get_twin_indices_rcpp(
    data = as.matrix(Xnorm),
    g = as.integer(g),
    v = as.integer(v),
    runs = as.integer(runs),
    seed = as.integer(seed)
  )
}