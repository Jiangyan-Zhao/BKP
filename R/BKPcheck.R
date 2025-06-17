#' Check if an object is a valid BKP model
#'
#' @description
#' This function verifies whether an object is a properly structured Beta Kernel Process (BKP) model.
#' It checks both the class inheritance and the presence of required components.
#'
#' @param object Any R object to be tested
#'
#' @return `TRUE` if the object meets all criteria for a valid BKP model, `FALSE` otherwise.
#'
#' @keywords internal
#' @export
#' @examples
#' set.seed(123)
#' n <- 100
#' x <- seq(-2, 2, length = n)
#' true_pi <- (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' m <- sample(100, n, replace = TRUE)
#' y <- rbinom(n, size = m, prob = true_pi)
#' df <- data.frame(x = x, y = y, m = m)
#' xx <- seq(-2, 2, length = 100) #new data points
#' model <- fit.bkp(df)
#' is.bkp(model)
is.bkp <- function(object) {
  # Verify class inheritance
  if (!inherits(x = object, what = "bkp")) return(FALSE)

  # List of required components
  required_components <- c("bestTheta", "minLoss", "kfun", "kernel.type", "loss.type", "y", "m",
                           "normalizedX", "originalX", "alpha0", "beta0")

  # Check presence of all components
  all(required_components %in% names(object))
}
