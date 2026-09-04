test_that("package exports the public QC interface", {
  exports <- getNamespaceExports("decontamSensitivity")
  expect_true(all(c(
    "run_decontam_qc",
    "run_decontam_threshold_sweep",
    "plot_flagged_taxa_reads",
    "plot_taxa_reads_before_after",
    "plot_prevalence_enrichment"
  ) %in% exports))
})

test_that("attach hook provides a concise startup message", {
  expect_message(
    decontamSensitivity:::.onAttach(
      system.file(package = "decontamSensitivity"),
      "decontamSensitivity"
    ),
    "run_decontam_qc"
  )
})
