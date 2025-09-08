# This test file validates the functionality of the `summary` S3 method for
# BKP objects. It checks for correct output structure, data types, and
# consistency with other functions like `predict`.

test_that("summary.BKP returns a well-structured summary object with correct values", {
  # Set a seed for reproducibility
  set.seed(123)

  # -------------------------------------------------------------------------
  # Setup: Create a real BKP model object using fit_BKP
  # -------------------------------------------------------------------------

  # Define true success probability function (1D)
  true_pi_fun <- function(x) {
    (1 + exp(-x^2) * cos(10 * (1 - exp(-x)) / (1 + exp(-x)))) / 2
  }

  n <- 50
  d <- 1
  Xbounds <- matrix(c(-2, 2), nrow = d)
  X <- tgp::lhs(n = n, rect = Xbounds)
  true_pi <- true_pi_fun(X)
  m <- sample(100, n, replace = TRUE)
  y <- rbinom(n, size = m, prob = true_pi)

  # Fit BKP model (this will be the object to test)
  model <- fit_BKP(X, y, m, Xbounds = Xbounds)

  # -------------------------------------------------------------------------
  # Test Cases
  # -------------------------------------------------------------------------

  # Case 1: Test basic summary properties
  summary_obj <- summary(model)

  # Check if the returned object is a list
  expect_type(summary_obj, "list")

  # Check if the object has the correct class for print method
  expect_s3_class(summary_obj, "summary_BKP")

  # Check for all expected named elements
  expected_names <- c("n_obs", "input_dim", "kernel", "theta_opt", "loss",
                      "loss_min", "prior", "r0", "p0", "post_mean", "post_var")
  expect_named(summary_obj, expected_names)

  # Case 2: Validate the values of the summary elements

  # Check number of observations
  expect_equal(summary_obj$n_obs, n)

  # Check input dimensionality
  expect_equal(summary_obj$input_dim, d)

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
  expect_equal(summary_obj$post_mean, as.vector(pred_result$mean))

  # The posterior variance from summary should be the same as from predict
  expect_equal(summary_obj$post_var, as.vector(pred_result$variance))

})

