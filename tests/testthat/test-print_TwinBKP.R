# This test file validates the functionality of the `print` S3 method for
# TwinBKP objects and their related output classes (summary, predict, simulate).
# It captures console output so that print methods do not pollute test output.

test_that("TwinBKP print methods produce expected output", {
  model <- make_twinbkp_model_1d()$model

  expect_output(
    print(model),
    "Twin Beta Kernel Process"
  )

  expect_output(
    print(predict(model)),
    "TwinBKP prediction"
  )

  expect_output(
    print(simulate(model, nsim = 2, seed = 1)),
    "TwinBKP Posterior Probability Simulations"
  )
})
