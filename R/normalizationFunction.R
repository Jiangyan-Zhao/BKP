#' @title Normalize a Matrix or Vector to the [0, 1] Range
#'
#' @description
#' This internal function scales the values of a numeric matrix or vector `X` to the [0, 1] range.
#' Each column is normalized independently. If a column contains constant values (i.e., its minimum
#' equals its maximum), all values in that column are set to 0.5 after normalization to
#' prevent division by zero and ensure a consistent representation.
#'
#' @param X A numeric matrix or vector to be normalized.
#'
#' @return A numeric matrix or vector with values scaled to the [0, 1] range.
#'   If the input was a vector, a numeric vector is returned. If the input was a matrix,
#'   a matrix is returned.
#'
#' @details
#' Normalization is performed using the formula:
#' \ifelse{html}{\out{<i>x<sub>norm</sub> = (x - min(X)) / (max(X) - min(X))</i>}}{\eqn{x_{\text{norm}} = \frac{x - \min(X)}{\max(X) - \min(X)}}}.
#' #' This process transforms the data to a common scale, which is often crucial for algorithms
#' sensitive to feature magnitudes, such as kernel methods. Special handling is implemented
#' for constant columns to ensure robust behavior.
#'
#' @keywords internal
#' This function is not exported for direct user use.
#'
normalize01 <- function(X) {
  # Step 1: Handle the case where X is a single column vector or a simple numeric vector.
  # If ncol(X) is NULL (for a vector) or 1 (for a single-column matrix), treat it as a single feature.
  if (base::is.null(base::ncol(X)) || base::ncol(X) == 1) {
    # Ensure X is treated as a matrix for consistent indexing, even if it was a simple vector.
    X_mat <- base::as.matrix(X)
    min_val <- base::min(X_mat, na.rm = TRUE) # Find the minimum value, ignoring NA.
    max_val <- base::max(X_mat, na.rm = TRUE) # Find the maximum value, ignoring NA.

    # Check if all values in the column are identical (min equals max).
    if (min_val == max_val) {
      # If constant, normalize all values to 0.5 to avoid division by zero and maintain a mid-range value.
      return(base::rep(0.5, base::nrow(X_mat))) # Return as a vector/column.
    } else {
      # Otherwise, perform standard min-max normalization to [0, 1].
      return((X_mat - min_val) / (max_val - min_val))
    }
  }

  # Step 2: Handle multi-column matrices.
  # Calculate the minimum and maximum values for each column.
  min_vals <- base::apply(X, 2, base::min, na.rm = TRUE)
  max_vals <- base::apply(X, 2, base::max, na.rm = TRUE)
  # Calculate the range (max - min) for each column.
  ranges <- max_vals - min_vals

  # Identify columns where the range is zero (i.e., all values are constant).
  const_cols <- base::which(ranges == 0)
  # To avoid division by zero during the `sweep` operation, temporarily set their range to 1.
  # These constant columns will be handled specifically afterwards.
  ranges[const_cols] <- 1

  # Perform normalization: (x - min) / range for each column.
  # `sweep(X, 2, min_vals, FUN = "-")` subtracts `min_vals` from each column of `X`.
  X_norm <- base::sweep(X, 2, min_vals, FUN = "-")
  # `sweep(X_norm, 2, ranges, FUN = "/")` divides each column of `X_norm` by its respective `ranges`.
  X_norm <- base::sweep(X_norm, 2, ranges, FUN = "/")

  # Step 3: Post-processing for constant columns.
  # For columns that were identified as constant, set all their normalized values to 0.5.
  if (base::length(const_cols) > 0) {
    X_norm[, const_cols] <- 0.5
  }

  return(X_norm) # Return the normalized matrix.
}
