test_that("quantile.TwinDKP returns DKP-style quantiles", { q <- quantile(make_twindkp_model_1d()$model, probs = c(0.025, 0.975)); expect_equal(dim(q), c(40, 3, 2)) })
