test_that("plot helpers return ggplot objects", {
  result <- run_toy_sweep()
  expect_s3_class(plot_threshold_sensitivity(result), "ggplot")
  expect_s3_class(plot_decontam_scores(result), "ggplot")
  expect_s3_class(plot_sample_read_retention(result), "ggplot")
})
