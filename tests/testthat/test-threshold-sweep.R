test_that("threshold sweep calculates read and feature retention", {
  result <- run_toy_sweep()
  expect_s3_class(result, "decontam_sensitivity")
  expect_equal(result$threshold_summary$contaminant_features, c(1, 2, 4))

  reads <- summarize_read_retention(result)
  at_01 <- reads[reads$threshold == 0.1 & reads$group == "biological", ]
  expect_equal(at_01$reads_before, 242)
  expect_equal(at_01$reads_retained, 232)

  features <- summarize_feature_retention(result)
  at_01_control <- features[
    features$threshold == 0.1 & features$group == "control", ]
  expect_equal(at_01_control$features_before, 4)
  expect_equal(at_01_control$features_retained, 3)
})

test_that("sample dropout is recorded", {
  result <- run_toy_sweep()
  at_04 <- result$sample_retention[result$sample_retention$threshold == 0.4, ]
  expect_true(all(at_04$zero_read[at_04$group == "control"]))
  expect_false(any(at_04$zero_read[at_04$group == "biological"]))

  summary <- summarize_sample_read_retention(result)
  expect_equal(nrow(summary), 6)
  at_depth <- summarize_samples_at_depth(result, depth = 50)
  expect_equal(nrow(at_depth), 3)
  expect_error(summarize_samples_at_depth(result, depth = 0), "positive")
})

test_that("orientation is inferred safely", {
  regular <- run_toy_sweep(FALSE)
  transposed <- run_toy_sweep(TRUE)
  expect_equal(regular$threshold_summary, transposed$threshold_summary)
})

test_that("missing decontam scores are not flagged", {
  result <- run_toy_sweep()
  asv6 <- result$feature_flags[result$feature_flags$feature_id == "ASV6", ]
  expect_false(any(asv6$contaminant))
})

test_that("invalid inputs fail explicitly", {
  x <- toy_sensitivity_data()
  bad_metadata <- x$metadata[-1, , drop = FALSE]
  expect_error(
    run_threshold_sweep(x$scores, x$counts, bad_metadata, "type", "control"),
    "sample IDs"
  )
  expect_error(
    run_threshold_sweep(x$scores, x$counts, x$metadata, "missing", "control"),
    "not present"
  )
  expect_error(
    run_threshold_sweep(x$scores, x$counts, x$metadata, "type", "control",
                        thresholds = 1.1),
    "from 0 to 1"
  )

  missing_group <- x$metadata
  missing_group$type[1] <- NA_character_
  expect_error(
    run_threshold_sweep(x$scores, x$counts, missing_group, "type", "control"),
    "missing values"
  )
})
