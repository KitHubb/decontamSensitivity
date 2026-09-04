#' Run a complete decontam threshold QC workflow
#'
#' One-call interface for a phyloseq object. It calculates prevalence scores,
#' evaluates all thresholds, filters at a prespecified threshold, summarizes
#' library size and control-versus-biological prevalence, and prepares the core
#' QC plots. The function does not select an optimal threshold.
#'
#' @param ps A phyloseq object.
#' @inheritParams run_decontam_threshold_sweep
#' @param selected_threshold Prespecified threshold used for taxonomic summaries
#'   and the filtered phyloseq object. It must occur in `thresholds`.
#' @param comparison_threshold_low Optional lower threshold for the newly flagged
#'   taxa table. If `NULL`, the nearest evaluated threshold below
#'   `selected_threshold` is used when available.
#' @param taxonomy Taxonomic rank used in summaries and plots.
#' @param top_n_flagged Maximum flagged taxa shown in the read-count plot.
#' @param top_n_composition Number of taxa shown individually in the before/after
#'   composition plot.
#' @param prevalence_pseudocount Correction used in prevalence enrichment ratios.
#' @param plot_all_thresholds Logical; when `TRUE`, also create named lists of
#'   flagged-taxa and before/after composition plots for every evaluated
#'   threshold. The selected-threshold plots are always returned separately.
#' @param progress Logical; display a `cli` progress bar. The default is `TRUE`
#'   in interactive sessions and `FALSE` during non-interactive execution such
#'   as R Markdown rendering and automated tests.
#' @return An object of class `decontam_qc` with `result`, `filtered_phyloseq`,
#'   `group_phyloseq`, `tables`, and `plots` components.
#' @export
#' @examples
#' \dontrun{
#' qc <- run_decontam_qc(
#'   ps,
#'   control_column = "control_status",
#'   control_label = "control",
#'   thresholds = seq(0.1, 0.9, 0.1),
#'   selected_threshold = 0.5,
#'   plot_all_thresholds = TRUE,
#'   progress = TRUE
#' )
#' qc$plots$sample_retention
#' qc$plots$taxa_reads_before_after
#' }
run_decontam_qc <- function(ps,
                            control_column,
                            control_label,
                            thresholds = seq(0.1, 0.5, 0.1),
                            selected_threshold,
                            comparison_threshold_low = NULL,
                            taxonomy = "Genus",
                            batch = NULL,
                            batch_combine = c("minimum", "product", "fisher"),
                            normalize = TRUE,
                            top_n_flagged = 20L,
                            top_n_composition = 15L,
                            prevalence_pseudocount = 0.5,
                            plot_all_thresholds = FALSE,
                            progress = interactive()) {
  thresholds <- .validate_thresholds(thresholds)
  if (!is.logical(plot_all_thresholds) ||
      length(plot_all_thresholds) != 1L || is.na(plot_all_thresholds)) {
    stop("`plot_all_thresholds` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (!is.numeric(selected_threshold) || length(selected_threshold) != 1L ||
      is.na(selected_threshold) || !any(abs(thresholds - selected_threshold) <
                                        sqrt(.Machine$double.eps))) {
    stop("`selected_threshold` must be one value included in `thresholds`.",
         call. = FALSE)
  }
  selected_threshold <- thresholds[
    which.min(abs(thresholds - selected_threshold))
  ]
  if (is.null(comparison_threshold_low)) {
    lower <- thresholds[thresholds < selected_threshold]
    comparison_threshold_low <- if (length(lower)) max(lower) else NULL
  }
  if (!is.null(comparison_threshold_low)) {
    if (!is.numeric(comparison_threshold_low) ||
        length(comparison_threshold_low) != 1L ||
        !any(abs(thresholds - comparison_threshold_low) <
             sqrt(.Machine$double.eps)) ||
        comparison_threshold_low >= selected_threshold) {
      stop("`comparison_threshold_low` must be an evaluated threshold below `selected_threshold`.",
           call. = FALSE)
    }
    comparison_threshold_low <- thresholds[
      which.min(abs(thresholds - comparison_threshold_low))
    ]
  }

  plot_steps <- if (isTRUE(plot_all_thresholds)) {
    2L * length(thresholds)
  } else {
    2L
  }
  progress_id <- NULL
  if (isTRUE(progress)) {
    progress_id <- cli::cli_progress_bar(
      name = "decontam QC",
      type = "tasks",
      total = 4L + plot_steps,
      clear = FALSE
    )
    on.exit({
      if (!is.null(progress_id)) {
        cli::cli_progress_done(id = progress_id)
      }
    }, add = TRUE)
  }
  update_progress <- function(status) {
    if (!is.null(progress_id)) {
      cli::cli_progress_update(
        id = progress_id, inc = 1L, status = status
      )
    }
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
  update_progress("decontam prevalence scores")
  filtered <- filter_phyloseq_at_threshold(
    ps, result, selected_threshold, prune_zero = TRUE
  )
  groups <- split_phyloseq_groups(
    ps, control_column, control_label, prune_zero = TRUE
  )
  update_progress("filtered phyloseq and sample groups")
  flagged_taxa <- summarize_flagged_taxa(
    result, selected_threshold, taxonomy
  )
  flag_rows <- result$feature_flags[
    result$feature_flags$threshold == selected_threshold &
      result$feature_flags$contaminant,
    ,
    drop = FALSE
  ]
  flagged_features <- result$feature_summary[
    match(flag_rows$feature_id, result$feature_summary$feature_id),
    ,
    drop = FALSE
  ]
  flagged_features <- flagged_features[
    order(-flagged_features$reads_biological, -flagged_features$reads_control),
    ,
    drop = FALSE
  ]
  newly_flagged <- NULL
  if (!is.null(comparison_threshold_low)) {
    newly_flagged <- compare_newly_flagged_taxa(
      result, comparison_threshold_low, selected_threshold, taxonomy
    )
  }
  tables <- list(
    threshold_summary = result$threshold_summary,
    read_retention = summarize_read_retention(result),
    feature_retention = summarize_feature_retention(result),
    sample_retention = result$sample_retention,
    sample_retention_summary = summarize_sample_read_retention(result),
    library_size_summary = summarize_library_sizes_by_group(
      ps, control_column, control_label
    ),
    prevalence_comparison = summarize_control_prevalence(
      ps, control_column, control_label, prevalence_pseudocount
    ),
    prevalence_enrichment = calculate_prevalence_enrichment(
      ps, control_column, control_label,
      prevalence_unit = "count", pseudocount = 0
    ),
    flagged_features = flagged_features,
    flagged_taxa = flagged_taxa,
    newly_flagged_taxa = newly_flagged,
    newly_flagged_intervals = summarize_newly_flagged_intervals(
      result, taxonomy
    )
  )
  update_progress("QC summary tables")
  threshold_names <- format(thresholds, trim = TRUE, scientific = FALSE)
  selected_name <- format(
    selected_threshold, trim = TRUE, scientific = FALSE
  )
  if (isTRUE(plot_all_thresholds)) {
    flagged_taxa_by_threshold <- stats::setNames(
      lapply(thresholds, function(threshold) {
        taxa_at_threshold <- summarize_flagged_taxa(
          result, threshold, taxonomy
        )
        plot_object <- if (nrow(taxa_at_threshold)) {
          plot_flagged_taxa_reads(
            result, threshold, taxonomy, top_n_flagged,
            measure = "mean_per_sample"
          )
        } else {
          NULL
        }
        update_progress(paste("flagged taxa at threshold", threshold))
        plot_object
      }),
      threshold_names
    )
    taxa_before_after_by_threshold <- stats::setNames(
      lapply(thresholds, function(threshold) {
        plot_object <- plot_taxa_reads_before_after(
          ps, result, threshold, control_column, control_label,
          taxonomy, top_n_composition
        )
        update_progress(paste("before/after taxa at threshold", threshold))
        plot_object
      }),
      threshold_names
    )
    flagged_taxa_plot <- flagged_taxa_by_threshold[[selected_name]]
    taxa_before_after_plot <-
      taxa_before_after_by_threshold[[selected_name]]
  } else {
    flagged_taxa_by_threshold <- NULL
    taxa_before_after_by_threshold <- NULL
    flagged_taxa_plot <- if (nrow(flagged_taxa)) {
      plot_flagged_taxa_reads(
        result, selected_threshold, taxonomy, top_n_flagged,
        measure = "mean_per_sample"
      )
    } else {
      NULL
    }
    update_progress(paste("flagged taxa at threshold", selected_threshold))
    taxa_before_after_plot <- plot_taxa_reads_before_after(
      ps, result, selected_threshold, control_column, control_label,
      taxonomy, top_n_composition
    )
    update_progress(paste("before/after taxa at threshold", selected_threshold))
  }
  plots <- list(
    score_distribution = plot_decontam_scores(result),
    threshold_sensitivity = plot_threshold_sensitivity(result),
    sample_retention = plot_sample_read_retention(result),
    prevalence_enrichment = plot_prevalence_enrichment(
      ps, control_column, control_label,
      prevalence_unit = "count", pseudocount = 0
    ),
    flagged_taxa_reads = flagged_taxa_plot,
    taxa_reads_before_after = taxa_before_after_plot,
    flagged_taxa_reads_by_threshold = flagged_taxa_by_threshold,
    taxa_reads_before_after_by_threshold = taxa_before_after_by_threshold
  )
  update_progress("core QC plots")
  if (!is.null(progress_id)) {
    cli::cli_progress_done(id = progress_id)
    progress_id <- NULL
  }
  structure(list(
    result = result,
    selected_threshold = selected_threshold,
    comparison_threshold_low = comparison_threshold_low,
    filtered_phyloseq = filtered,
    group_phyloseq = groups,
    tables = tables,
    plots = plots,
    call = match.call()
  ), class = "decontam_qc")
}

#' @export
print.decontam_qc <- function(x, ...) {
  cat("Complete decontam threshold QC\n")
  cat("  selected threshold: ", x$selected_threshold, "\n", sep = "")
  cat("  retained features: ", phyloseq::ntaxa(x$filtered_phyloseq), "\n", sep = "")
  cat("  retained samples: ", phyloseq::nsamples(x$filtered_phyloseq), "\n", sep = "")
  if (!is.null(x$comparison_threshold_low)) {
    cat("  newly flagged comparison: ", x$comparison_threshold_low,
        " -> ", x$selected_threshold, "\n", sep = "")
  }
  invisible(x)
}
