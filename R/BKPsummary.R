#' @title Summarize a BKP Model
#'
#' @description
#' This function provides a summary of a fitted Beta Kernel Process (BKP) model.
#' Currently, it acts as a wrapper around the `print.bkp` method, displaying
#' the same concise overview of the model's key characteristics and fitting results.
#'
#' @param object An object of class "bkp", typically returned by the `fit.bkp()` function.
#' @param ... Additional arguments for compatibility with the generic `summary` method
#'   (currently not used by this specific method).
#'
#' @return This function is called for its side effect of printing a summary to the console.
#'   It returns `NULL` invisibly.
#'
#' @details
#' The `summary.bkp` method leverages the existing `print.bkp` function to provide
#' essential information about the fitted BKP model. This includes the number of
#' observations, covariate dimensions, kernel type, optimized parameters, and the
#' minimum loss value. In future versions, this function might be expanded to
#' provide more detailed statistical summaries or diagnostic information.
#'
#' @seealso
#' \code{\link{print.bkp}} for the basic print method.
#' \code{\link{fit.bkp}} to fit a BKP model.
#'
#' @export
#' @method summary bkp
#'
#' @examples
#' ### 1D
#' set.seed(123)
#' n <- 100
#' x <- seq(-2, 2, length = n)
#' true_pi <- (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#' df <- data.frame(x = x, y = y, m = m)
#' xx <- seq(-2, 2, length = 100) #new data points
#' model <- fit.bkp(df)
#' summary(model)
summary.bkp <- function(object, ...){
  # Validate the input object: Ensure it is a valid BKP model.
  # This check is implicitly handled by print.bkp, but can be added here for directness if preferred.
  # if (!inherits(object, "bkp")) {
  #   stop("The object in question is not of class \"bkp\". Please provide a model fitted by 'fit.bkp()'.")
  # }

  # Call the print.bkp method to display the model summary.
  # This provides a consistent output format for both print and summary.
  print.bkp(object, ...)

  invisible(NULL) # Return NULL invisibly, as the function's primary purpose is printing.
}
