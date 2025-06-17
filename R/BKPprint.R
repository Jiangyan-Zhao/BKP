#' @title Print a BKP Model Summary
#'
#' @description
#' This function provides a concise summary of a fitted Beta Kernel Process (BKP) model
#' object when it's printed to the console. It displays key information such as
#' the number of observations, the dimensions of covariates, the kernel type used,
#' the optimized kernel parameters (`bestTheta`), and the minimum Brier score achieved.
#'
#' @param x An object of class "bkp", typically returned by the `fit.bkp()` function.
#' @param ... Additional arguments for compatibility with the generic `print` method
#'   (currently not used by this specific method).
#'
#' @return This function is called for its side effect of printing a summary to the console.
#'   It returns `NULL` invisibly.
#'
#' @details
#' When a `bkp` object is passed to `print()` (or simply typed at the R console),
#' this method will be invoked to present a structured overview of the model's key characteristics
#' and fitting results. It helps in quickly assessing the fitted model.
#'
#' @seealso
#' \code{\link{fit.bkp}} to fit a BKP model.
#' \code{\link{summary.bkp}} for a more detailed summary if available.
#'
#' @export
#' @method print bkp
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
#' print(model)
print.bkp <- function(x, ...) {
  bkp.model <- x # Assign the input 'x' (the bkp object) to a more descriptive variable.

  # Validate the input object: Ensure it is a valid BKP model.
  # Using 'inherits' for robust class checking.
  if (!inherits(bkp.model, "bkp")) {
    base::stop("The object in question is not of class \"bkp\". Please provide a model fitted by 'fit.bkp()'.")
  }

  # Print a formatted header for the BKP model summary.
  base::cat("---\n")
  base::cat("## BKP Model Summary\n")
  base::cat("---\n\n")

  # Print key summary statistics of the BKP model.
  # 'nrow(bkp.model$y)' gives the number of observations (n).
  base::cat("Number of Observations (n): ", base::nrow(bkp.model$y), "\n", sep = "")
  # 'ncol(bkp.model$originalX)' gives the number of covariate dimensions (d).
  base::cat("X Input Dimensions (d): ", base::ncol(bkp.model$originalX), "\n", sep = "")
  # Display the type of kernel used in the model.
  base::cat("Kernel Type: ", bkp.model$kernel.type, "\n", sep = "")
  # Display the type of loss used in the model.
  base::cat("Loss Type: ", bkp.model$loss.type, "\n", sep = "")
  # Show the optimized kernel parameters (bestTheta), rounded for readability.
  base::cat("Optimized Best Theta: ", base::paste(base::round(bkp.model$bestTheta, 4), collapse = ", "), "\n", sep = "")
  # Display the minimum loss achieved during model fitting, rounded for readability.
  base::cat("Corresponding Minimum Loss : ", base::round(bkp.model$minLoss, 4), "\n", sep = "")

  invisible(NULL) # Return NULL invisibly, as the function's primary purpose is printing.
}
