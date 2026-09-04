#' Extract read-retention results
#'
#' @param result A result from [run_threshold_sweep()].
#' @return A data frame with one row per threshold and sample group.
#' @export
summarize_read_retention <- function(result) {
  if (!inherits(result, "decontam_sensitivity")) {
    stop("`result` must be returned by `run_threshold_sweep()`.", call. = FALSE)
  }
  result$read_retention
}

#' Extract feature-retention results
#'
#' Group-specific feature denominators include only features observed in that
#' group before filtering. This avoids counting control-only features as
#' biological features, and vice versa.
#'
#' @inheritParams summarize_read_retention
#' @return A data frame with one row per threshold and sample group.
#' @export
summarize_feature_retention <- function(result) {
  if (!inherits(result, "decontam_sensitivity")) {
    stop("`result` must be returned by `run_threshold_sweep()`.", call. = FALSE)
  }
  result$feature_retention
}

#' Summarize sample-level read retention
#'
#' Converts the per-sample results into threshold-by-group summaries suitable
#' for reporting alongside [plot_sample_read_retention()].
#'
#' @inheritParams summarize_read_retention
#' @return A data frame with sample count, zero-read count, and the minimum and
#'   median percentage of reads retained for each threshold and sample group.
#' @export
summarize_sample_read_retention <- function(result) {
  if (!inherits(result, "decontam_sensitivity")) {
    stop("`result` must be returned by `run_threshold_sweep()`.", call. = FALSE)
  }
  do.call(rbind, lapply(result$thresholds, function(threshold) {
    do.call(rbind, lapply(c("biological", "control"), function(group_name) {
      rows <- result$sample_retention$threshold == threshold &
        result$sample_retention$group == group_name
      values <- result$sample_retention$reads_retained_pct[rows]
      data.frame(
        threshold = threshold,
        group = group_name,
        samples = sum(rows),
        zero_read_samples = sum(result$sample_retention$zero_read[rows]),
        minimum_reads_retained_pct = if (all(is.na(values))) NA_real_ else min(values, na.rm = TRUE),
        median_reads_retained_pct = if (all(is.na(values))) NA_real_ else stats::median(values, na.rm = TRUE),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }))
  }))
}

#' Summarize samples remaining above a read-depth cutoff
#'
#' This is an optional downstream-eligibility diagnostic, not a decontam
#' threshold-selection criterion.
#'
#' @inheritParams summarize_read_retention
#' @param depth Positive read-depth cutoff.
#' @param group One of `"biological"` or `"control"`.
#' @return A threshold-level data frame.
#' @export
summarize_samples_at_depth <- function(result,
                                       depth,
                                       group = c("biological", "control")) {
  if (!inherits(result, "decontam_sensitivity")) {
    stop("`result` must be returned by `run_threshold_sweep()`.", call. = FALSE)
  }
  if (!is.numeric(depth) || length(depth) != 1L || is.na(depth) ||
      !is.finite(depth) || depth <= 0) {
    stop("`depth` must be one positive finite number.", call. = FALSE)
  }
  group <- match.arg(group)
  do.call(rbind, lapply(result$thresholds, function(threshold) {
    rows <- result$sample_retention$threshold == threshold &
      result$sample_retention$group == group
    reads <- result$sample_retention$reads_retained[rows]
    data.frame(
      threshold = threshold,
      group = group,
      samples = length(reads),
      samples_available_at_depth = sum(reads >= depth),
      samples_below_depth = sum(reads < depth),
      minimum_retained_reads = min(reads),
      median_retained_reads = stats::median(reads),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
}

#' Summarize taxonomy flagged at one threshold
#'
#' @inheritParams summarize_read_retention
#' @param threshold A threshold included in the sweep.
#' @param taxonomy Taxonomy column to aggregate, such as `"Genus"` or
#'   `"Species"`.
#' @return A taxonomy summary ordered by biological-sample read abundance.
#' @export
summarize_flagged_taxa <- function(result, threshold, taxonomy = "Genus") {
  i <- .threshold_index(result, threshold)
  taxa <- .taxonomy_for_result(result, taxonomy)
  flagged <- result$feature_flags$threshold == result$thresholds[i] &
    result$feature_flags$contaminant
  ids <- result$feature_flags$feature_id[flagged]
  .aggregate_taxa(result, ids, taxa)
}

#' Compare taxa newly flagged between two thresholds
#'
#' Reports features classified as contaminants at `threshold_high` but not at
#' `threshold_low`. The result makes biological-signal loss visible without
#' claiming that either threshold is optimal.
#'
#' @inheritParams summarize_flagged_taxa
#' @param threshold_low Lower threshold included in the sweep.
#' @param threshold_high Higher threshold included in the sweep.
#' @return A taxonomy summary for newly flagged features.
#' @export
compare_newly_flagged_taxa <- function(result,
                                       threshold_low,
                                       threshold_high,
                                       taxonomy = "Genus") {
  low <- .threshold_index(result, threshold_low)
  high <- .threshold_index(result, threshold_high)
  if (result$thresholds[low] >= result$thresholds[high]) {
    stop("`threshold_low` must be less than `threshold_high`.", call. = FALSE)
  }
  taxa <- .taxonomy_for_result(result, taxonomy)
  low_flags <- result$feature_flags$feature_id[
    result$feature_flags$threshold == result$thresholds[low] &
      result$feature_flags$contaminant
  ]
  high_flags <- result$feature_flags$feature_id[
    result$feature_flags$threshold == result$thresholds[high] &
      result$feature_flags$contaminant
  ]
  .aggregate_taxa(result, setdiff(high_flags, low_flags), taxa)
}

#' Summarize taxa newly flagged across adjacent thresholds
#'
#' @inheritParams summarize_flagged_taxa
#' @return A data frame combining every adjacent threshold comparison.
#' @export
summarize_newly_flagged_intervals <- function(result, taxonomy = "Genus") {
  if (!inherits(result, "decontam_sensitivity")) {
    stop("`result` must be returned by `run_threshold_sweep()`.", call. = FALSE)
  }
  .taxonomy_for_result(result, taxonomy)
  empty <- data.frame(
    threshold_low = numeric(), threshold_high = numeric(),
    taxon = character(), newly_flagged_features = integer(),
    reads_biological = numeric(),
    relative_abundance_biological_pct = numeric(),
    reads_control = numeric(), relative_abundance_control_pct = numeric(),
    stringsAsFactors = FALSE
  )
  if (length(result$thresholds) < 2L) return(empty)
  pieces <- lapply(seq.int(2L, length(result$thresholds)), function(i) {
    out <- compare_newly_flagged_taxa(
      result, result$thresholds[i - 1L], result$thresholds[i], taxonomy
    )
    if (!nrow(out)) return(NULL)
    out$threshold_low <- result$thresholds[i - 1L]
    out$threshold_high <- result$thresholds[i]
    out[, c(
      "threshold_low", "threshold_high", "taxon",
      "newly_flagged_features", "reads_biological",
      "relative_abundance_biological_pct", "reads_control",
      "relative_abundance_control_pct"
    )]
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) return(empty)
  rownames_out <- do.call(rbind, pieces)
  rownames(rownames_out) <- NULL
  rownames_out
}

.aggregate_taxa <- function(result, feature_ids, taxa) {
  empty <- data.frame(
    taxon = character(), newly_flagged_features = integer(),
    reads_biological = numeric(), relative_abundance_biological_pct = numeric(),
    reads_control = numeric(), relative_abundance_control_pct = numeric(),
    stringsAsFactors = FALSE
  )
  if (!length(feature_ids)) return(empty)

  fs <- result$feature_summary
  rows <- match(feature_ids, fs$feature_id)
  x <- data.frame(
    taxon = taxa[rows], feature_id = feature_ids,
    reads_biological = fs$reads_biological[rows],
    reads_control = fs$reads_control[rows],
    stringsAsFactors = FALSE
  )
  pieces <- split(x, x$taxon, drop = TRUE)
  out <- do.call(rbind, lapply(pieces, function(z) {
    data.frame(
      taxon = z$taxon[1L],
      newly_flagged_features = nrow(z),
      reads_biological = sum(z$reads_biological),
      relative_abundance_biological_pct = .safe_percent(
        sum(z$reads_biological), sum(result$counts[, !result$is_control, drop = FALSE])
      ),
      reads_control = sum(z$reads_control),
      relative_abundance_control_pct = .safe_percent(
        sum(z$reads_control), sum(result$counts[, result$is_control, drop = FALSE])
      ),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out[order(-out$reads_biological, -out$newly_flagged_features, out$taxon), , drop = FALSE]
}
