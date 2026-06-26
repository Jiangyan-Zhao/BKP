test_that("fitted.TwinDKP returns fitted probabilities", { fit <- make_twindkp_model_1d()$model; expect_equal(fitted(fit), fit$prob) })
