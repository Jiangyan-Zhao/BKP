#' @title Calculate Weighted Euclidean Distance
#'
#' @description
#' This internal function computes the weighted Euclidean distance between a single vector `x1`
#' and another vector or a matrix `x2`. The distance calculation is influenced by `theta`,
#' a vector of weights that scales the contribution of each dimension.
#'
#' @param x1 A numeric vector representing a single data point. Its length defines the number of dimensions.
#' @param x2 A numeric vector or a matrix.
#'   If `x2` is a vector, it is treated as a single data point to which `x1` is compared.
#'   If `x2` is a matrix, each row is considered a separate data point, and the distance
#'   between `x1` and each row of `x2` is calculated.
#' @param theta A numeric vector of the same length as `x1` (and the number of columns of `x2` if `x2` is a matrix).
#'   These values act as weights for each dimension. The distance is computed using the formula:
#'   \ifelse{html}{\out{ <i>d</i> = &radic; [ &Sigma;<sub>i=1</sub><sup>D</sup> ((x1<sub>i</sub> - x2<sub>i</sub>) / &theta;<sub>i</sub>)<sup>2</sup> ] }}{\deqn{
#'     d = \sqrt{\sum_{i=1}^{D} \left( \frac{x1_i - x2_i}{\theta_i} \right)^2}
#'   }}
#'   where \eqn{D} is the number of dimensions.
#'
#' @return A numeric value if `x2` is a vector, representing the single weighted Euclidean distance
#'   between `x1` and `x2`. If `x2` is a matrix, it returns a numeric vector where each element
#'   is the weighted Euclidean distance between `x1` and the corresponding row of `x2`.
#'
#' @details
#' This function is a core component for kernel-based methods where the distance metric
#' needs to be adaptable by different scaling parameters for each dimension.
#' It's typically used internally within kernel functions to compute similarities between data points.
#'
#' @keywords internal
#' # This function is not exported for direct user use.
distFun <- function(x1, x2, theta) {
  # Check if x2 is a single vector (i.e., not a matrix, as indicated by NULL for nrow).
  # This distinguishes between calculating distance to a single point or to multiple points.
  if (base::is.null(base::nrow(x2))) {
    # If x2 is a vector, calculate the weighted Euclidean distance to this single point.
    # Step 1: Calculate squared differences for each dimension, scaled by theta.
    diff_squared <- ((x1 - x2) / theta)^2
    # Step 2: Sum these squared differences and take the square root.
    return(base::sqrt(base::sum(diff_squared)))
  } else {
    # If x2 is a matrix, calculate the weighted Euclidean distance from x1 to each row of x2.
    # `apply()` iterates over each row (MARGIN = 1) of x2.
    distances <- base::apply(x2, MARGIN = 1, FUN = function(row_vec) {
      # For each row, perform the same weighted Euclidean distance calculation.
      diff_squared <- ((x1 - row_vec) / theta)^2
      return(base::sqrt(base::sum(diff_squared)))
    })
    return(distances) # Return a vector of distances.
  }
}
