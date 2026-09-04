.require_phyloseq <- function() {
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package `phyloseq` is required for this function.", call. = FALSE)
  }
}

.phyloseq_group_info <- function(ps, control_column, control_label) {
  .require_phyloseq()
  if (!inherits(ps, "phyloseq")) {
    stop("`ps` must be a phyloseq object.", call. = FALSE)
  }
  .assert_scalar_character(control_column, "control_column")
  metadata <- as.data.frame(phyloseq::sample_data(ps))
  if (!control_column %in% names(metadata)) {
    stop("Column `", control_column, "` is not present in sample data.",
         call. = FALSE)
  }
  if (anyNA(metadata[[control_column]])) {
    stop("`control_column` contains missing values; classify every sample explicitly.",
         call. = FALSE)
  }
  if (!length(control_label) || anyNA(control_label)) {
    stop("`control_label` must contain at least one non-missing value.",
         call. = FALSE)
  }
  is_control <- metadata[[control_column]] %in% control_label
  if (!any(is_control) || all(is_control)) {
    stop("Both biological samples and negative controls are required.",
         call. = FALSE)
  }
  names(is_control) <- rownames(metadata)
  list(metadata = metadata, is_control = is_control)
}

.phyloseq_count_matrix <- function(ps) {
  counts <- methods::as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) counts <- t(counts)
  counts
}

#' Calculate abundance and prevalence for every feature
#'
#' This is the package version of the original project helper with the same
#' output column names. Abundance is the total raw-read count and prevalence is
#' the number of samples in which the feature is detected.
#'
#' @param ps A phyloseq object.
#' @return A data frame with taxa as row names and columns `Abundance` and
#'   `Prevalence`.
#' @export
get_abundance_prevalence <- function(ps) {
  .require_phyloseq()
  if (!inherits(ps, "phyloseq")) {
    stop("`ps` must be a phyloseq object.", call. = FALSE)
  }
  counts <- .phyloseq_count_matrix(ps)
  data.frame(
    Abundance = as.numeric(phyloseq::taxa_sums(ps)[rownames(counts)]),
    Prevalence = rowSums(counts > 0),
    row.names = rownames(counts),
    check.names = FALSE
  )
}

#' Calculate sample-versus-control prevalence enrichment
#'
#' Reproduces the `kA`, `odds.sample`, `odds.control`, and `summarized`
#' calculations from the original analysis in one call. Despite the legacy
#' column name, `odds.sample` is a prevalence ratio, not a statistical odds
#' ratio.
#'
#' @inheritParams split_phyloseq_groups
#' @param prevalence_unit Use the original prevalence counts (`"count"`) or
#'   prevalence proportions (`"proportion"`). Proportions are preferable when
#'   the two groups contain different numbers of samples.
#' @param pseudocount Non-negative correction added to the ratio inputs. Use
#'   zero with `prevalence_unit = "count"` to reproduce the original code
#'   exactly.
#' @return A feature-level data frame containing the original project column
#'   names plus explicit count and proportion columns.
#' @export
calculate_prevalence_enrichment <- function(ps,
                                             control_column,
                                             control_label,
                                             prevalence_unit = c("count", "proportion"),
                                             pseudocount = 0) {
  prevalence_unit <- match.arg(prevalence_unit)
  if (!is.numeric(pseudocount) || length(pseudocount) != 1L ||
      is.na(pseudocount) || !is.finite(pseudocount) || pseudocount < 0) {
    stop("`pseudocount` must be one finite non-negative number.", call. = FALSE)
  }
  groups <- split_phyloseq_groups(
    ps, control_column, control_label, prune_zero = FALSE
  )
  control <- get_abundance_prevalence(groups$control)
  sample <- get_abundance_prevalence(groups$biological)
  feature_ids <- phyloseq::taxa_names(ps)
  control <- control[feature_ids, , drop = FALSE]
  sample <- sample[feature_ids, , drop = FALSE]

  control_prop <- control$Prevalence / phyloseq::nsamples(groups$control)
  sample_prop <- sample$Prevalence / phyloseq::nsamples(groups$biological)
  if (prevalence_unit == "count") {
    ratio_control <- control$Prevalence
    ratio_sample <- sample$Prevalence
  } else {
    ratio_control <- control_prop
    ratio_sample <- sample_prop
  }

  data.frame(
    feature_id = feature_ids,
    Abundance.Control = control$Abundance,
    Prevalence.Control = ratio_control,
    Abundance.Sample = sample$Abundance,
    Prevalence.Sample = ratio_sample,
    odds.sample = (ratio_sample + pseudocount) /
      (ratio_control + pseudocount),
    odds.control = (ratio_control + pseudocount) /
      (ratio_sample + pseudocount),
    summarized = sample$Abundance + control$Abundance,
    Prevalence.Count.Control = control$Prevalence,
    Prevalence.Count.Sample = sample$Prevalence,
    Prevalence.Proportion.Control = control_prop,
    Prevalence.Proportion.Sample = sample_prop,
    stringsAsFactors = FALSE,
    row.names = feature_ids,
    check.names = FALSE
  )
}

#' Split a phyloseq object into biological and control samples
#'
#' Replaces repeated `subset_samples()` and `prune_taxa()` blocks with one
#' explicit, reusable operation.
#'
#' @param ps A phyloseq object.
#' @param control_column Sample-data column identifying negative controls.
#' @param control_label Value or values identifying negative controls.
#' @param prune_zero Logical; remove features with zero reads within each group.
#' @return A named list containing `biological` and `control` phyloseq objects.
#' @export
split_phyloseq_groups <- function(ps,
                                  control_column,
                                  control_label,
                                  prune_zero = TRUE) {
  group_info <- .phyloseq_group_info(ps, control_column, control_label)
  biological <- phyloseq::prune_samples(!group_info$is_control, ps)
  control <- phyloseq::prune_samples(group_info$is_control, ps)
  if (isTRUE(prune_zero)) {
    biological <- phyloseq::prune_taxa(
      phyloseq::taxa_sums(biological) > 0, biological
    )
    control <- phyloseq::prune_taxa(
      phyloseq::taxa_sums(control) > 0, control
    )
  }
  list(biological = biological, control = control)
}

#' Summarize library sizes by sample group
#'
#' @inheritParams split_phyloseq_groups
#' @return A data frame containing group size and the standard library-size
#'   summary statistics.
#' @export
summarize_library_sizes_by_group <- function(ps,
                                             control_column,
                                             control_label) {
  groups <- split_phyloseq_groups(
    ps, control_column, control_label, prune_zero = FALSE
  )
  do.call(rbind, lapply(names(groups), function(group_name) {
    sizes <- as.numeric(phyloseq::sample_sums(groups[[group_name]]))
    data.frame(
      group = group_name,
      samples = length(sizes),
      minimum = min(sizes),
      first_quartile = unname(stats::quantile(sizes, 0.25)),
      median = stats::median(sizes),
      mean = mean(sizes),
      third_quartile = unname(stats::quantile(sizes, 0.75)),
      maximum = max(sizes),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
}

#' Summarize feature abundance and prevalence by group
#'
#' Uses read totals for abundance and the proportion of samples containing each
#' feature for prevalence. Corrected prevalence ratios use a small pseudocount
#' so that features absent from one group remain finite and plottable.
#'
#' @inheritParams split_phyloseq_groups
#' @param pseudocount Non-negative correction used for prevalence proportions.
#'   With the default 0.5, the corrected proportion is
#'   `(observed + 0.5) / (samples + 1)`.
#' @return A feature-level data frame aligned to `taxa_names(ps)`.
#' @export
summarize_control_prevalence <- function(ps,
                                         control_column,
                                         control_label,
                                         pseudocount = 0.5) {
  group_info <- .phyloseq_group_info(ps, control_column, control_label)
  if (!is.numeric(pseudocount) || length(pseudocount) != 1L ||
      is.na(pseudocount) || !is.finite(pseudocount) || pseudocount < 0) {
    stop("`pseudocount` must be one finite non-negative number.", call. = FALSE)
  }
  counts <- .phyloseq_count_matrix(ps)
  is_control <- group_info$is_control[colnames(counts)]
  n_biological <- sum(!is_control)
  n_control <- sum(is_control)
  prevalence_biological <- rowSums(counts[, !is_control, drop = FALSE] > 0)
  prevalence_control <- rowSums(counts[, is_control, drop = FALSE] > 0)
  prevalence_pct_biological <- 100 * prevalence_biological / n_biological
  prevalence_pct_control <- 100 * prevalence_control / n_control

  if (pseudocount == 0) {
    corrected_biological <- prevalence_biological / n_biological
    corrected_control <- prevalence_control / n_control
  } else {
    corrected_biological <-
      (prevalence_biological + pseudocount) /
      (n_biological + 2 * pseudocount)
    corrected_control <-
      (prevalence_control + pseudocount) /
      (n_control + 2 * pseudocount)
  }

  data.frame(
    feature_id = rownames(counts),
    abundance_biological = rowSums(counts[, !is_control, drop = FALSE]),
    abundance_control = rowSums(counts[, is_control, drop = FALSE]),
    abundance_total = rowSums(counts),
    prevalence_biological = prevalence_biological,
    prevalence_control = prevalence_control,
    prevalence_pct_biological = prevalence_pct_biological,
    prevalence_pct_control = prevalence_pct_control,
    prevalence_ratio_biological = corrected_biological / corrected_control,
    prevalence_ratio_control = corrected_control / corrected_biological,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Filter a phyloseq object at one evaluated threshold
#'
#' @param result A result from [run_threshold_sweep()].
#' @param threshold A threshold included in `result`.
#' @inheritParams split_phyloseq_groups
#' @return A filtered phyloseq object.
#' @export
filter_phyloseq_at_threshold <- function(ps,
                                         result,
                                         threshold,
                                         prune_zero = TRUE) {
  .require_phyloseq()
  i <- .threshold_index(result, threshold)
  counts <- .phyloseq_count_matrix(ps)
  if (!setequal(rownames(counts), rownames(result$counts)) ||
      !setequal(colnames(counts), colnames(result$counts))) {
    stop("`ps` and `result` must contain the same feature and sample IDs.",
         call. = FALSE)
  }
  flags <- result$feature_flags
  remove_ids <- flags$feature_id[
    flags$threshold == result$thresholds[i] & flags$contaminant
  ]
  keep <- !phyloseq::taxa_names(ps) %in% remove_ids
  filtered <- phyloseq::prune_taxa(keep, ps)
  if (isTRUE(prune_zero)) {
    filtered <- phyloseq::prune_taxa(
      phyloseq::taxa_sums(filtered) > 0, filtered
    )
  }
  filtered
}
