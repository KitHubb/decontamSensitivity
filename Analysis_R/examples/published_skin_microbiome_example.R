# Published skin microbiome workflow
#
# Replace the path below with the phyloseq object used for the Scientific
# Reports analysis. The original study data are not redistributed here.

library(decontamSensitivity)

# ps <- readRDS("path/to/published_skin_phyloseq.rds")
# result <- run_decontam_threshold_sweep(
#   ps = ps,
#   control_column = "T_C",
#   control_label = "Ctrl",
#   thresholds = seq(0.1, 0.5, 0.1),
#   batch = "SequencingRun"
# )
#
# result$threshold_summary
# plot_threshold_sensitivity(result)
# plot_decontam_scores(result)
# summarize_flagged_taxa(result, threshold = 0.1, taxonomy = "Genus")
