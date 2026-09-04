test_that("newly flagged taxa are summarized with biological abundance", {
  result <- run_toy_sweep()
  out <- compare_newly_flagged_taxa(result, 0.1, 0.2, "Genus")
  expect_equal(nrow(out), 1)
  expect_equal(out$taxon, "Staphylococcus")
  expect_equal(out$newly_flagged_features, 1)
  expect_equal(out$reads_biological, 10)
  expect_equal(out$reads_control, 10)
})

test_that("all flagged taxa can be summarized", {
  result <- run_toy_sweep()
  out <- summarize_flagged_taxa(result, 0.2, "Genus")
  expect_equal(out$newly_flagged_features[out$taxon == "Staphylococcus"], 2)
  expect_error(compare_newly_flagged_taxa(result, 0.4, 0.1), "must be less")
  expect_error(summarize_flagged_taxa(result, 0.3), "not included")

  intervals <- summarize_newly_flagged_intervals(result, "Genus")
  expect_equal(sort(unique(intervals$threshold_high)), c(0.2, 0.4))
})
