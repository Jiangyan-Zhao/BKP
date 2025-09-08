# Test file for the kernel_matrix function

test_that("kernel_matrix handles input validation correctly", {
  # Test with different column numbers for X and Xprime
  expect_error(kernel_matrix(matrix(1:4, 2), matrix(1:6, 2)),
               "'X' and 'Xprime' must have the same number of columns \\(input dimensions\\).")

  # Test with an invalid kernel name
  expect_error(kernel_matrix(matrix(1:4, 2), kernel = "invalid_kernel"),
               "should be one of \"gaussian\", \"matern52\", \"matern32\"", fixed = TRUE)

  # Test with anisotropic=FALSE and theta is a vector
  expect_error(kernel_matrix(matrix(1:4, 2), anisotropic = FALSE, theta = c(0.5, 0.6)),
               "For anisotropic=FALSE, 'theta' must be a scalar.")

  # Test with anisotropic=TRUE and theta has wrong length
  expect_error(kernel_matrix(matrix(1:6, 3), anisotropic = TRUE, theta = c(0.5, 0.6, 0.7)),
               "For anisotropic=TRUE, 'theta' must be scalar or of length equal to ncol\\(X\\).")
})


