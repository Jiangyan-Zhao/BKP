#' @title Predict New Data Points using a BKP Model
#'
#' @description
#' This function generates predictions for new data points using a previously fitted
#' Beta Kernel Process (BKP) model.
#'
#' @param object A fitted BKP model object (typically obtained from `fit.bkp()`).
#' @param Xnew New data points for prediction. Should be a matrix or data frame.
#' @param CI.size The significance level for prediction intervals (default: 0.05).
#'
#' @return A data frame with prediction results:
#' \item{X}{Original input new data points}
#' \item{piHat}{Predicted probabilities of success}
#' \item{lower}{Lower bound of confidence interval}
#' \item{upper}{Upper bound of confidence interval}
#'
#' @details
#' The prediction process involves:
#' \enumerate{
#'   \item Model validation
#'   \item Data normalization
#'   \item Kernel matrix calculation
#'   \item Posterior parameter calculation
#'   \item Prediction and CI computation
#' }
#'
#' @seealso
#' \code{\link{fit.bkp}} to fit a BKP model.
#' \code{\link{plot.bkp}} to visualize predictions.
#'
#' @method predict bkp
#' @export
#' @importFrom stats qbeta
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
#' head(predict(model,xx))

predict.bkp <- function(object, Xnew, CI.size = 0.05) {
  bkp.model <- object # Assign the input 'object' to a more descriptive variable.

  # Step 1: Validate the input object.
  # Ensure that the provided 'object' is a properly fitted BKP model.
  if (!is.bkp(bkp.model)) {
    stop("The object in question is not of class \"bkp\". Please provide a model fitted by 'fit.bkp()'.")
  }

  # Step 2: Prepare new data points for prediction.
  # Convert Xnew to a matrix and store the original version before normalization.
  OriginalXnew <- Xnew <- as.matrix(Xnew)
  # Normalize Xnew using the same scaling as the training data.
  # (Assumes 'normalize01' is an internal helper function that scales X to [0, 1]).
  Xnew <- normalize01(Xnew)

  # Step 3: Compute the Kernel Matrix (K).
  # K represents the similarity between each new data point in 'Xnew' and
  # every normalized data point from the original training data (bkp.model$normalizedX).
  # The 'kfun' (kernel function) defined in the BKP model is used here.
  K <- apply(Xnew, MARGIN = 1, FUN = function(x1_row_current) {
    bkp.model$kfun(x1_row_current, bkp.model$normalizedX, bkp.model$bestTheta)
  })
  # Transpose K to align dimensions for matrix multiplication (N_new x N_original).
  K <- t(K)

  # Step 4: Calculate Posterior Parameters (Alpha and Beta).
  # These are the shape parameters for the posterior Beta distributions
  # for each new data point, based on the prior parameters (alpha0, beta0)
  # and weighted observations (y, m-y) through the kernel matrix.
  Alpha <- bkp.model$alpha0 + K %*% bkp.model$y
  Beta <- bkp.model$beta0 + K %*% (bkp.model$m - bkp.model$y)

  # Step 5: Calculate Predicted Probabilities and Confidence Intervals.
  # piHat is the mean of the posterior Beta distribution (E[p] = alpha / (alpha + beta)).
  piHat <- Alpha / (Alpha + Beta)
  # piLower and piUpper are calculated using quantiles of the posterior Beta distribution.
  # qbeta(p, shape1, shape2) gives the p-th quantile of the Beta distribution.
  # CI.size/2 for the lower bound, 1 - CI.size/2 for the upper bound.
  piLower <- stats::qbeta(CI.size / 2, Alpha, Beta)
  piUpper <- stats::qbeta(1 - CI.size / 2, Alpha, Beta)

  # Return the results as a data frame.
  # The 'X' column contains the original (unnormalized) new data points.
  return(data.frame(X = OriginalXnew, piHat = piHat, lower = piLower, upper = piUpper))
}
