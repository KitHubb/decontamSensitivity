#' Run a complete decontam threshold QC workflow
#'
#' One-call interface for a phyloseq object. It calculates prevalence scores,
#' evaluates every supplied threshold, filters the object at every threshold,
#' summarizes QC metrics, and prepares threshold-specific plots. The function
#' does not select or compare a preferred threshold.
#'
#' @param ps A phyloseq object.
#' @inheritParams run_decontam_threshold_sweep
#' @param taxonomy Taxonomic rank used in summaries and plots.
#' @param top_n_flagged Maximum flagged taxa shown in each read-count plot.
#' @param top_n_composition Number of taxa shown individually in each
#'   before/after composition plot.
#' @param prevalence_pseudocount Correction used in prevalence enrichment ratios.
#' @param group_colors Named colors for `total`, `biological`, and `control`.
#' @param taxa_colors Optional named colors for taxa. `NULL` applies the
#'   automatic RColorBrewer palette rule.
#' @param progress Logical; display one updating `cli` progress bar. The default
#'   is `TRUE` in interactive sessions and `FALSE` during non-interactive
#'   execution such as R Markdown rendering and automated tests.
#' @return An object of class `decontam_qc` containing results, a named list of
#'   filtered phyloseq objects, summary tables, and threshold-specific plots.
#' @export
#' @examples
#' \dontrun{
#' qc <- run_decontam_qc(
#'   ps,
#'   control_column = "control_status",
#'   control_label = "control",
#'   thresholds = seq(0.1, 0.9, 0.1),
#'   progress = TRUE
#' )
#' qc$plots$sample_retention
#' qc$plots$taxa_reads_before_after_by_threshold[["0.5"]]
#' }
run_decontam_qc <- function(ps,
                            control_column,
                            control_label,
                            thresholds = seq(0.1, 0.5, 0.1),
                            taxonomy = "Genus",
                            batch = NULL,
                            batch_combine = c("minimum", "product", "fisher"),
                            normalize = TRUE,
                            top_n_flagged = 20L,
                            top_n_composition = 15L,
                            prevalence_pseudocount = 0.5,
                            group_colors = c(
                              total = "black",
                              biological = "steelblue",
                              control = "tomato"
                            ),
                            taxa_colors = NULL,
                            progress = interactive()) {
  thresholds <- .validate_thresholds(thresholds)
  group_colors <- .validate_group_colors(group_colors)
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  progress_id <- NULL
  if (isTRUE(progress)) {
    progress_id <- cli::cli_progress_bar(
      name = "decontam QC",
      total = 4L + 3L * length(thresholds),
      format = "{cli::pb_name} {cli::pb_bar} {cli::pb_percent}",
      clear = FALSE
    )
    on.exit({
      if (!is.null(progress_id)) cli::cli_progress_done(id = progress_id)
    }, add = TRUE)
  }
  update_progress <- function() {
    if (!is.null(progress_id)) cli::cli_progress_update(id = progress_id)
    invisible(NULL)
  }

  result <- run_decontam_threshold_sweep(
    ps = ps,
    control_column = control_column,
    control_label = control_label,
    thresholds = thresholds,
    batch = batch,
    batch_combine = match.arg(batch_combine),
    normalize = normalize
  )
  update_progress()
  threshold_names <- format(thresholds, trim = TRUE, scientific = FALSE)
  filtered_by_threshold <- stats::setNames(
    lapply(thresholds, function(threshold) {
      out <- filter_phyloseq_at_threshold(ps, result, threshold, prune_zero = TRUE)
      update_progress()
      out
    }),
    threshold_names
  )
  groups <- split_phyloseq_groups(ps, control_column, control_label, prune_zero = TRUE)
  update_progress()

  flagged_feature_parts <- lapply(thresholds, function(threshold) {
    flags <- result$feature_flags[
      result$feature_flags$threshold == threshold & result$feature_flags$contaminant,
      , drop = FALSE
    ]
    out <- result$feature_summary[
      match(flags$feature_id, result$feature_summary$feature_id),
      , drop = FALSE
    ]
    out$threshold <- threshold
    out[, c("threshold", setdiff(names(out), "threshold")), drop = FALSE]
  })
  flagged_features <- do.call(rbind, flagged_feature_parts)
  rownames(flagged_features) <- NULL

  flagged_taxa_parts <- lapply(thresholds, function(threshold) {
    out <- summarize_flagged_taxa(result, threshold, taxonomy)
    out$threshold <- threshold
    out[, c("threshold", setdiff(names(out), "threshold")), drop = FALSE]
  })
  flagged_taxa <- do.call(rbind, flagged_taxa_parts)
  rownames(flagged_taxa) <- NULL

  tables <- list(
    threshold_summary = result$threshold_summary,
    read_retention = summarize_read_retention(result),
    feature_retention = summarize_feature_retention(result),
    sample_retention = result$sample_retention,
    sample_retention_summary = summarize_sample_read_retention(result),
    library_size_summary = summarize_library_sizes_by_group(ps, control_column, control_label),
    prevalence_comparison = summarize_control_prevalence(
      ps, control_column, control_label, prevalence_pseudocount
    ),
    prevalence_enrichment = calculate_prevalence_enrichment(
      ps, control_column, control_label, prevalence_unit = "count", pseudocount = 0
    ),
    flagged_features = flagged_features,
    flagged_taxa = flagged_taxa
  )
  update_progress()

  flagged_taxa_by_threshold <- stats::setNames(
    lapply(thresholds, function(threshold) {
      taxa_at_threshold <- summarize_flagged_taxa(result, threshold, taxonomy)
      out <- if (nrow(taxa_at_threshold)) {
        plot_flagged_taxa_reads(
          result, threshold, taxonomy, top_n_flagged,
          measure = "mean_per_sample", group_colors = group_colors
        )
      } else NULL
      update_progress()
      out
    }),
    threshold_names
  )
  taxa_before_after_by_threshold <- stats::setNames(
    lapply(thresholds, function(threshold) {
      out <- plot_taxa_reads_before_after(
        ps, result, threshold, control_column, control_label,
        taxonomy, top_n_composition, taxa_colors = taxa_colors
      )
      update_progress()
      out
    }),
    threshold_names
  )
  plots <- list(
    score_distribution = plot_decontam_scores(result),
    threshold_sensitivity = plot_threshold_sensitivity(
      result, group_colors = group_colors
    ),
    sample_retention = plot_sample_read_retention(
      result, group_colors = group_colors
    ),
    prevalence_enrichment = plot_prevalence_enrichment(
      ps, control_column, control_label, prevalence_unit = "count", pseudocount = 0
    ),
    flagged_taxa_reads_by_threshold = flagged_taxa_by_threshold,
    taxa_reads_before_after_by_threshold = taxa_before_after_by_threshold
  )
  update_progress()
  if (!is.null(progress_id)) {
    cli::cli_progress_done(id = progress_id)
    progress_id <- NULL
  }

  structure(list(
    result = result,
    thresholds = thresholds,
    filtered_phyloseq_by_threshold = filtered_by_threshold,
    group_phyloseq = groups,
    tables = tables,
    plots = plots,
    call = match.call()
  ), class = "decontam_qc")
}

#' @export
print.decontam_qc <- function(x, ...) {
  cat("Complete decontam threshold QC\n")
  cat("  thresholds: ", paste(x$thresholds, collapse = ", "), "\n", sep = "")
  cat("  filtered phyloseq objects: ", length(x$filtered_phyloseq_by_threshold),
      "\n", sep = "")
  invisible(x)
}
