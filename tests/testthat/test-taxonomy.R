test_that("all flagged taxa can be summarized", {
  result <- run_toy_sweep()
  out <- summarize_flagged_taxa(result, 0.2, "Genus")
  expect_equal(out$flagged_features[out$taxon == "Staphylococcus"], 2)
  expect_error(summarize_flagged_taxa(result, 0.3), "not included")
})
