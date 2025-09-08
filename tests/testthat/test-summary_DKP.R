# This test file validates the functionality of the `summary` S3 method for
# DKP objects. It checks for correct output structure, data types, and
# consistency with other functions like `predict`.

test_that("summary.DKP returns a well-structured summary object with correct values", {
  # Set a seed for reproducibility
  set.seed(123)

  # -------------------------------------------------------------------------
  # Setup: Create a real DKP model object using fit_DKP
  # -------------------------------------------------------------------------

  # Define true class probability function (3-class)
  true_pi_fun <- function(X) {
    p1 <- 1/(1+exp(-3*X))
    p2 <- (1 + exp(-X^2) * cos(10 * (1 - exp(-X)) / (1 + exp(-X)))) / 2
    return(matrix(c(p1/2, p2/2, 1 - (p1+p2)/2), nrow = length(p1)))
  }

  n <- 30
  d <- 1
  q <- 3
  Xbounds <- matrix(c(-2, 2), nrow = d)
  X <- tgp::lhs(n = n, rect = Xbounds)
  true_pi <- true_pi_fun(X)
  m <- sample(150, n, replace = TRUE)

  # Generate multinomial responses
  Y <- t(sapply(1:n, function(i) rmultinom(1, size = m[i], prob = true_pi[i, ])))

  # Fit DKP model (this will be the object to test)
  model <- fit_DKP(X, Y, Xbounds = Xbounds)

  # -------------------------------------------------------------------------
  # Test Cases
  # -------------------------------------------------------------------------

  # Case 1: Test basic summary properties
  summary_obj <- summary(model)

  # Check if the returned object is a list
  expect_type(summary_obj, "list")

  # Check if the object has the correct class for print method
  expect_s3_class(summary_obj, "summary_DKP")

  # Check for all expected named elements, adjusted for the current function's output
  # NOTE: The test now expects 'n_class' instead of 'output_dim'
  expected_names <- c("n_obs", "input_dim", "n_class", "kernel", "theta_opt", "loss",
                      "loss_min", "prior", "r0", "p0", "post_mean", "post_var")
  expect_named(summary_obj, expected_names)

  # Case 2: Validate the values of the summary elements

  # Check number of observations
  expect_equal(summary_obj$n_obs, n)

  # Check input and output dimensionality
  expect_equal(summary_obj$input_dim, d)
  # NOTE: The test now checks 'n_class' which is the correct name from the function's output
  expect_equal(summary_obj$n_class, q)

  # Dynamically check kernel and loss type from the fitted model object
  expect_equal(summary_obj$kernel, model$kernel)
  expect_equal(summary_obj$loss, model$loss)

  # Check if theta_opt is a numeric vector
  expect_type(summary_obj$theta_opt, "double")
  expect_true(all(summary_obj$theta_opt >= 0))

  # Check if loss_min is a single numeric value
  expect_type(summary_obj$loss_min, "double")
  expect_length(summary_obj$loss_min, 1)

  # Dynamically check prior parameters from the fitted model object
  expect_equal(summary_obj$prior, model$prior)
  expect_equal(summary_obj$r0, model$r0)
  expect_equal(summary_obj$p0, model$p0)

  # Case 3: Verify consistency with `predict` function

  # Get predictions on the training data using the predict method
  pred_result <- predict(model)

  # The posterior mean from summary should be the same as from predict
  expect_equal(summary_obj$post_mean, as.matrix(pred_result$mean))

  # The posterior variance from summary should be the same as from predict
  expect_equal(summary_obj$post_var, as.matrix(pred_result$variance))

})


