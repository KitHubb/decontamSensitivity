test_that("plot helpers return ggplot objects", {
  result <- run_toy_sweep()
  expect_s3_class(plot_threshold_sensitivity(result), "ggplot")
  expect_s3_class(plot_decontam_scores(result), "ggplot")
  expect_s3_class(plot_sample_read_retention(result), "ggplot")
})

test_that("standard group colors and taxon palettes are applied", {
  result <- run_toy_sweep()
  p <- plot_threshold_sensitivity(result)
  color_scale <- p$scales$get_scales("colour")
  expect_equal(
    unname(color_scale$palette(3)),
    c("black", "steelblue", "tomato")
  )
  expect_length(decontamSensitivity:::.taxa_palette(9), 9)
  expect_length(decontamSensitivity:::.taxa_palette(12), 12)
  expect_equal(
    decontamSensitivity:::.taxa_palette(21),
    c(
      RColorBrewer::brewer.pal(12, "Paired"),
      RColorBrewer::brewer.pal(9, "Set1")
    )
  )
})
