make_toy_ps <- function() {
  x <- toy_sensitivity_data()
  phyloseq::phyloseq(
    phyloseq::otu_table(x$counts, taxa_are_rows = TRUE),
    phyloseq::sample_data(x$metadata),
    phyloseq::tax_table(as.matrix(x$taxonomy))
  )
}

test_that("phyloseq convenience wrapper runs one prevalence model", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("decontam")
  ps <- make_toy_ps()
  result <- suppressWarnings(run_decontam_threshold_sweep(
    ps, control_column = "type", control_label = "control",
    thresholds = c(0.1, 0.2)
  ))
  expect_s3_class(result, "decontam_sensitivity")
  expect_equal(result$thresholds, c(0.1, 0.2))
})

test_that("phyloseq groups and library sizes are summarized", {
  skip_if_not_installed("phyloseq")
  ps <- make_toy_ps()
  groups <- split_phyloseq_groups(ps, "type", "control")
  expect_equal(phyloseq::nsamples(groups$biological), 2)
  expect_equal(phyloseq::nsamples(groups$control), 2)

  libraries <- summarize_library_sizes_by_group(
    ps, "type", "control"
  )
  expect_equal(libraries$group, c("biological", "control"))
  expect_equal(libraries$samples, c(2, 2))
})

test_that("control prevalence summary uses comparable proportions", {
  skip_if_not_installed("phyloseq")
  ps <- make_toy_ps()
  summary <- summarize_control_prevalence(
    ps, "type", "control"
  )
  expect_equal(summary$feature_id, phyloseq::taxa_names(ps))
  expect_true(all(is.finite(summary$prevalence_ratio_biological)))
  expect_true(all(is.finite(summary$prevalence_ratio_control)))
  expect_s3_class(
    plot_prevalence_enrichment(ps, "type", "control"),
    "ggplot"
  )
})

test_that("original prevalence enrichment block is available as functions", {
  skip_if_not_installed("phyloseq")
  ps <- make_toy_ps()
  one_group <- split_phyloseq_groups(ps, "type", "control")
  helper <- get_abundance_prevalence(one_group$control)
  expect_named(helper, c("Abundance", "Prevalence"))

  enrichment <- calculate_prevalence_enrichment(
    ps, "type", "control", prevalence_unit = "count"
  )
  expect_true(all(c(
    "Abundance.Control", "Prevalence.Control", "Abundance.Sample",
    "Prevalence.Sample", "odds.sample", "odds.control", "summarized"
  ) %in% names(enrichment)))
  expect_equal(
    enrichment$odds.sample,
    enrichment$Prevalence.Sample / enrichment$Prevalence.Control
  )
  expect_s3_class(
    plot_prevalence_enrichment(ps, "type", "control"), "ggplot"
  )
})

test_that("phyloseq helpers give identical results for both OTU orientations", {
  skip_if_not_installed("phyloseq")
  ps_taxa_rows <- make_toy_ps()
  counts <- methods::as(phyloseq::otu_table(ps_taxa_rows), "matrix")
  ps_samples_rows <- phyloseq::phyloseq(
    phyloseq::otu_table(t(counts), taxa_are_rows = FALSE),
    phyloseq::sample_data(ps_taxa_rows),
    phyloseq::tax_table(ps_taxa_rows)
  )

  expect_false(phyloseq::taxa_are_rows(ps_samples_rows))
  expect_equal(
    get_abundance_prevalence(ps_taxa_rows),
    get_abundance_prevalence(ps_samples_rows)
  )
  expect_equal(
    calculate_prevalence_enrichment(
      ps_taxa_rows, "type", "control", "count", 0
    ),
    calculate_prevalence_enrichment(
      ps_samples_rows, "type", "control", "count", 0
    )
  )
  expect_equal(
    summarize_control_prevalence(ps_taxa_rows, "type", "control"),
    summarize_control_prevalence(ps_samples_rows, "type", "control")
  )
})

test_that("phyloseq-specific plots and filtering are reusable", {
  skip_if_not_installed("phyloseq")
  x <- toy_sensitivity_data()
  ps <- make_toy_ps()
  result <- run_threshold_sweep(
    x$scores,
    x$counts,
    x$metadata,
    "type",
    "control",
    thresholds = c(0.1, 0.5),
    taxonomy = x$taxonomy
  )
  filtered <- filter_phyloseq_at_threshold(ps, result, 0.5)
  expect_lt(phyloseq::ntaxa(filtered), phyloseq::ntaxa(ps))
  expect_s3_class(plot_sample_read_retention(result), "ggplot")
  expect_s3_class(plot_flagged_taxa_reads(result, 0.5), "ggplot")
  expect_s3_class(
    plot_taxa_reads_before_after(
      ps, result, 0.5, "type", "control"
    ),
    "ggplot"
  )
})

test_that("one-call QC returns tables, plots, and filtered data", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("decontam")
  ps <- make_toy_ps()
  qc <- suppressWarnings(run_decontam_qc(
    ps,
    control_column = "type",
    control_label = "control",
    thresholds = c(0.1, 0.2, 0.5)
  ))
  expect_s3_class(qc, "decontam_qc")
  expect_s3_class(qc$result, "decontam_sensitivity")
  expect_equal(qc$thresholds, c(0.1, 0.2, 0.5))
  expect_named(qc$filtered_phyloseq_by_threshold, c("0.1", "0.2", "0.5"))
  expect_true(all(vapply(qc$filtered_phyloseq_by_threshold,
                         methods::is, logical(1), class2 = "phyloseq")))
  expect_named(qc$plots, c(
    "score_distribution", "threshold_sensitivity", "sample_retention",
    "prevalence_enrichment", "flagged_taxa_reads_by_threshold",
    "taxa_reads_before_after_by_threshold"
  ))
  expect_named(qc$tables, c(
    "threshold_summary", "read_retention", "feature_retention",
    "sample_retention", "sample_retention_summary", "library_size_summary",
    "prevalence_comparison", "prevalence_enrichment", "flagged_features", "flagged_taxa"
  ))
  expect_true(all(qc$tables$flagged_features$threshold %in% qc$thresholds))
  expect_false(any(c("selected_threshold", "comparison_threshold_low",
                     "plot_all_thresholds") %in% names(formals(run_decontam_qc))))
})

test_that("one-call QC creates taxa plots for every threshold", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("decontam")
  ps <- make_toy_ps()
  qc <- suppressWarnings(run_decontam_qc(
    ps,
    control_column = "type",
    control_label = "control",
    thresholds = c(0.1, 0.2, 0.5)
  ))

  expect_named(
    qc$plots$flagged_taxa_reads_by_threshold,
    c("0.1", "0.2", "0.5")
  )
  expect_named(
    qc$plots$taxa_reads_before_after_by_threshold,
    c("0.1", "0.2", "0.5")
  )
  expect_true(all(vapply(
    qc$plots$taxa_reads_before_after_by_threshold,
    inherits, logical(1), what = "ggplot"
  )))
})

test_that("one-call QC validates the progress option", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("decontam")
  ps <- make_toy_ps()
  expect_error(
    run_decontam_qc(
      ps,
      control_column = "type",
      control_label = "control",
      thresholds = c(0.1, 0.5),
      progress = NA
    ),
    "`progress` must be `TRUE` or `FALSE`",
    fixed = TRUE
  )
})
