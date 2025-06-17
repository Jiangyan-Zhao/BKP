#' @title Calculate Confusion Matrix Proportions for a BKP Model
#'
#' @description
#' This function computes the proportions of true positives (TP), true negatives (TN),
#' false positives (FP), and false negatives (FN) based on the predictions
#' from a fitted Bayesian Kernel Probit (BKP) model. It is specifically designed for binary
#' classification tasks where the outcome `y` is binary (0 or 1).
#'
#' @param bkp.model An object of class "bkp", typically returned by the `fit.bkp()` function.
#' @param threshold A numeric value between 0 and 1. This is the cutoff point
#'   used to convert the predicted probabilities from the BKP model into binary
#'   classifications (0 or 1). Defaults to 0.5.
#'
#' @return A list containing the proportions of:
#' \itemize{
#'   \item \code{TP}: True Positives (proportion of actual 1s correctly predicted as 1).
#'   \item \code{TN}: True Negatives (proportion of actual 0s correctly predicted as 0).
#'   \item \code{FP}: False Positives (proportion of actual 0s incorrectly predicted as 1).
#'   \item \code{FN}: False Negatives (proportion of actual 1s incorrectly predicted as 0).
#' }
#' The sum of TP, TN, FP, and FN will always be 1.
#'
#' @details
#' The function first validates the input `bkp.model` and ensures that the
#' observed outcome `y` in the model is indeed binary (composed only of 0s and 1s).
#' It then uses the `predict.bkp` method to obtain predicted probabilities for the
#' original covariate data `X` from the model. These probabilities are converted
#' into binary predictions (0 or 1) by comparing them against the specified `threshold`.
#' Finally, it calculates the proportions of TP, TN, FP, and FN relative to the
#' total number of observations.
#'
#' @export
#'
#' @examples
#' ### 1D
#' set.seed(123)
#' n <- 100
#' x <- seq(-2, 2, length = n)
#' true_pi <- (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
#' m <- rep(1,n)
#' y <- rbinom(n, size = m, prob = true_pi)
#' df <- data.frame(x = x, y = y, m = m)
#' xx <- seq(-2, 2, length = 100) #new data points
#' model <- fit.bkp(df)
#' confusionCal.bkp(model)

confusionCal.bkp <- function(bkp.model, threshold = 0.5) {
  # Step 1: Input Validation - Check if the provided object is a valid "bkp" model.
  # This ensures that the function operates on the correct type of object.
  if (!is.bkp(bkp.model)) {
    base::stop("The object in question is not of class \"bkp\". Please provide a model fitted by 'fit.bkp()'.")
  }

  # Extract the observed outcome variable 'y' from the BKP model.
  y <- bkp.model$y

  # Step 2: Input Validation - Check if the outcome variable 'y' is binary.
  # This function is specifically for binary classification, so 'y' must contain only 0s and 1s.
  # 'unique(y)' gets distinct values, 'sort()' orders them, and 'all() == c(0, 1)' checks if they are exactly 0 and 1.
  if (!base::all(base::sort(base::unique(y)) == c(0, 1))) {
    base::stop("The outcome 'y' in the bkp.model is not a binary response (i.e., not composed of only 0s and 1s).")
  }

  # Extract the original covariate data 'X' from the BKP model.
  x <- bkp.model$originalX
  # Step 3: Predict probabilities using the fitted model on the original covariate data.
  # This uses the 'predict.bkp' method to get the predicted probabilities (piHat) for the training data.
  pred_probs <- predict.bkp(bkp.model, x)$piHat

  # Step 4: Convert predicted probabilities to binary classifications.
  # If the predicted probability is greater than the 'threshold', classify as 1; otherwise, classify as 0.
  # 'as.numeric()' converts the logical result (TRUE/FALSE) into 1/0.
  pred_y <- base::as.numeric(pred_probs > threshold)

  # Get the total number of observations.
  n <- base::length(y)

  # Step 5: Calculate the proportions for each component of the confusion matrix.
  # True Positives (TP): Cases where actual y is 1 AND predicted y is 1.
  TP <- base::sum(y == 1 & pred_y == 1) / n
  # True Negatives (TN): Cases where actual y is 0 AND predicted y is 0.
  TN <- base::sum(y == 0 & pred_y == 0) / n
  # False Positives (FP): Cases where actual y is 0 BUT predicted y is 1 (Type I error).
  FP <- base::sum(y == 0 & pred_y == 1) / n
  # False Negatives (FN): Cases where actual y is 1 BUT predicted y is 0 (Type II error).
  FN <- base::sum(y == 1 & pred_y == 0) / n

  # Step 6: Return the calculated proportions as a named list.
  return(base::list(TP = TP, TN = TN, FP = FP, FN = FN))
}
