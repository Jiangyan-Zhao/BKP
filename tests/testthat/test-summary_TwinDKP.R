test_that("summary.TwinDKP returns a summary object", { expect_s3_class(summary(make_twindkp_model_1d()$model), "summary_TwinDKP") })
