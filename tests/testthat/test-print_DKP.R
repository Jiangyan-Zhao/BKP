# This test file validates the functionality of the `print` S3 method for
# DKP objects and their related output classes (summary, predict, simulate).
# It captures console output so that print methods do not pollute test output.

test_that("print.DKP methods run without errors and produce expected output", {
  fit <- make_dkp_model_1d()
  model <- fit$model

  expect_output(
    print(model),
    "Dirichlet Kernel Process"
  )

  expect_output(
    print(summary(model)),
    "Dirichlet Kernel Process"
  )

  expect_output(
    print(predict(model)),
    "Prediction results"
  )

  expect_output(
    print(simulate(model)),
    "Simulation results"
  )
})


test_that("print.DKP handles new-data, multi-class and high-dimensional branches", {
  fit <- make_dkp_model_3d_classification()
  model <- fit$model
  Xnew <- fit$X[1:5, , drop = FALSE]

  expect_output(
    print(model),
    "Dirichlet Kernel Process"
  )

  expect_output(
    print(summary(model)),
    "Dirichlet Kernel Process"
  )

  expect_output(
    print(predict(model, Xnew = Xnew)),
    "Prediction results"
  )

  expect_output(
    print(simulate(model, Xnew = Xnew, nsim = 5)),
    "Simulation results"
  )
})
