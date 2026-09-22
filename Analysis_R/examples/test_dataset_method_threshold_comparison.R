# Frequency, prevalence, and combined threshold comparison on the test dataset
#
# This example follows the package tutorial and uses the small dataset from
# tests/testthat/helper-toy-data.R. It is intended to verify method and
# threshold behavior, not to select a threshold for a biological study.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package `devtools` is required to load the development package.")
}
if (!requireNamespace("phyloseq", quietly = TRUE)) {
  stop("Package `phyloseq` is required for this example.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required for this example.")
}

devtools::load_all(".", quiet = TRUE)
source(file.path("tests", "testthat", "helper-toy-data.R"))

toy <- toy_sensitivity_data()
ps_test <- phyloseq::phyloseq(
  phyloseq::otu_table(toy$counts, taxa_are_rows = TRUE),
  phyloseq::sample_data(toy$metadata),
  phyloseq::tax_table(as.matrix(toy$taxonomy))
)

thresholds <- c(0.01, 0.05, seq(0.1, 0.9, by = 0.1))
methods <- c("frequency", "prevalence", "combined")
method_labels <- c(
  frequency = "Frequency",
  prevalence = "Prevalence",
  combined = "Combined"
)

run_method <- function(method) {
  args <- list(
    ps = ps_test,
    control_column = "type",
    control_label = "control",
    thresholds = thresholds,
    method = method
  )
  if (method %in% c("frequency", "combined")) {
    args$concentration_column <- "DNA_concentration"
  }
  suppressWarnings(do.call(run_decontam_threshold_sweep, args))
}

results <- stats::setNames(lapply(methods, run_method), methods)

threshold_summary <- do.call(rbind, lapply(methods, function(method) {
  out <- results[[method]]$threshold_summary
  out$method <- unname(method_labels[method])
  out[, c("method", setdiff(names(out), "method")), drop = FALSE]
}))
rownames(threshold_summary) <- NULL

feature_flags <- do.call(rbind, lapply(methods, function(method) {
  out <- results[[method]]$feature_flags
  out$method <- unname(method_labels[method])
  out$genus <- toy$taxonomy[out$feature_id, "Genus"]
  out[, c("method", "threshold", "feature_id", "genus", "contaminant")]
}))
rownames(feature_flags) <- NULL

method_pairs <- combn(unname(method_labels[methods]), 2L, simplify = FALSE)
agreement_parts <- list()
part_index <- 0L
for (threshold in thresholds) {
  flags_at_threshold <- feature_flags[
    feature_flags$threshold == threshold,
    c("method", "feature_id", "contaminant")
  ]
  for (pair in method_pairs) {
    a <- flags_at_threshold[flags_at_threshold$method == pair[1L], ]
    b <- flags_at_threshold[flags_at_threshold$method == pair[2L], ]
    b <- b[match(a$feature_id, b$feature_id), ]
    both_flagged <- sum(a$contaminant & b$contaminant)
    either_flagged <- sum(a$contaminant | b$contaminant)
    part_index <- part_index + 1L
    agreement_parts[[part_index]] <- data.frame(
      threshold = threshold,
      method_1 = pair[1L],
      method_2 = pair[2L],
      features_compared = nrow(a),
      both_flagged = both_flagged,
      either_flagged = either_flagged,
      agreement_pct = 100 * mean(a$contaminant == b$contaminant),
      flagged_jaccard_pct = if (either_flagged == 0L) {
        NA_real_
      } else {
        100 * both_flagged / either_flagged
      }
    )
  }
}
method_agreement <- do.call(rbind, agreement_parts)

metric_columns <- c(
  contaminant_features = "Contaminant features (n)",
  biological_reads_retained_pct = "Biological reads retained (%)",
  control_reads_retained_pct = "Control reads retained (%)",
  biological_features_retained_pct = "Biological features retained (%)",
  control_features_retained_pct = "Control features retained (%)"
)
metric_parts <- lapply(names(metric_columns), function(column) {
  data.frame(
    method = threshold_summary$method,
    threshold = threshold_summary$threshold,
    metric = unname(metric_columns[column]),
    value = threshold_summary[[column]],
    stringsAsFactors = FALSE
  )
})
comparison_long <- do.call(rbind, metric_parts)
comparison_long$method <- factor(
  comparison_long$method,
  levels = unname(method_labels[methods])
)
comparison_long$metric <- factor(
  comparison_long$metric,
  levels = unname(metric_columns)
)

comparison_plot <- ggplot2::ggplot(
  comparison_long,
  ggplot2::aes(x = threshold, y = value, color = method, group = method)
) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 2) +
  ggplot2::scale_x_continuous(breaks = thresholds) +
  ggplot2::scale_color_manual(values = c(
    Frequency = "#0072B2",
    Prevalence = "#D55E00",
    Combined = "#009E73"
  )) +
  ggplot2::labs(
    title = "Test dataset: decontam method and threshold comparison",
    subtitle = "Frequency, prevalence, and combined methods across 11 thresholds",
    x = "decontam threshold",
    y = NULL,
    color = "Method"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

analysis_results_dir <- Sys.getenv(
  "DECONTAM_ANALYSIS_RESULTS",
  unset = file.path(
    "..", "decontamSensitivity_paper", "R_analysis", "Results"
  )
)
output_dir <- file.path(
  analysis_results_dir,
  "test_method_comparison"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

utils::write.csv(
  threshold_summary,
  file.path(output_dir, "threshold_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  feature_flags,
  file.path(output_dir, "feature_flags.csv"),
  row.names = FALSE
)
utils::write.csv(
  method_agreement,
  file.path(output_dir, "method_agreement.csv"),
  row.names = FALSE
)
ggplot2::ggsave(
  file.path(output_dir, "method_threshold_comparison.png"),
  comparison_plot,
  width = 10,
  height = 8,
  dpi = 300
)
saveRDS(
  list(
    thresholds = thresholds,
    results = results,
    threshold_summary = threshold_summary,
    feature_flags = feature_flags,
    method_agreement = method_agreement,
    plot = comparison_plot
  ),
  file.path(output_dir, "method_threshold_comparison.rds")
)

print(threshold_summary)
print(method_agreement)
if (interactive()) print(comparison_plot)
