#' Run a threshold sensitivity sweep
#'
#' Applies multiple classification thresholds to feature scores already
#' produced by [decontam::isContaminant()]. For each threshold, the function
#' quantifies read retention, feature retention, sample dropout, and feature
#' classifications in biological samples and negative controls.
#'
#' @param decontam_result A data frame containing one row per feature and a
#'   score column, normally the detailed output from
#'   [decontam::isContaminant()]. Feature IDs must be row names unless
#'   `feature_id_column` is supplied.
#' @param count_table Numeric feature-by-sample matrix. A sample-by-feature
#'   matrix is also accepted when `features_are_rows = FALSE`, or can be
#'   inferred from metadata IDs.
#' @param metadata Data frame with one row per sample. Sample IDs must be row
#'   names unless `sample_id_column` is supplied.
#' @param control_column Metadata column identifying negative controls.
#' @param control_label Value or values in `control_column` that identify
#'   negative controls.
#' @param thresholds Numeric thresholds between zero and one. A feature is
#'   classified as a contaminant when its score is strictly less than the
#'   threshold, matching `decontam`.
#' @param taxonomy Optional data frame containing feature taxonomy. Feature IDs
#'   must be row names unless `taxonomy_id_column` is supplied.
#' @param score_column Name of the contaminant score column.
#' @param feature_id_column Optional feature-ID column in `decontam_result`.
#' @param sample_id_column Optional sample-ID column in `metadata`.
#' @param taxonomy_id_column Optional feature-ID column in `taxonomy`.
#' @param features_are_rows Logical indicating count-table orientation. If
#'   `NULL`, orientation is inferred from sample IDs.
#'
#' @return An object of class `decontam_sensitivity` containing threshold,
#'   feature, read, and sample summaries plus aligned source data.
#' @export
#' @examples
#' \dontrun{
#' result <- run_threshold_sweep(
#'   decontam_result = decontam_scores,
#'   count_table = feature_counts,
#'   metadata = sample_metadata,
#'   control_column = "sample_class",
#'   control_label = "control",
#'   thresholds = seq(0.1, 0.5, 0.1),
#'   taxonomy = feature_taxonomy
#' )
#' result
#' }
run_threshold_sweep <- function(decontam_result,
                                count_table,
                                metadata,
                                control_column,
                                control_label,
                                thresholds = seq(0.1, 0.5, 0.1),
                                taxonomy = NULL,
                                score_column = "p",
                                feature_id_column = NULL,
                                sample_id_column = NULL,
                                taxonomy_id_column = NULL,
                                features_are_rows = NULL) {
  .assert_scalar_character(control_column, "control_column")
  .assert_scalar_character(score_column, "score_column")
  thresholds <- .validate_thresholds(thresholds)

  metadata <- .ids_from_data_frame(metadata, sample_id_column, "metadata")
  if (!control_column %in% names(metadata)) {
    stop("Column `", control_column, "` is not present in `metadata`.",
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

  counts <- as.matrix(count_table)
  if (!is.numeric(counts) || length(dim(counts)) != 2L || anyNA(counts) ||
      any(!is.finite(counts)) || any(counts < 0)) {
    stop("`count_table` must be a finite, non-negative numeric matrix.",
         call. = FALSE)
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
      anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
    stop("`count_table` must have unique row and column names.", call. = FALSE)
  }

  sample_ids <- rownames(metadata)
  row_match <- setequal(rownames(counts), sample_ids)
  column_match <- setequal(colnames(counts), sample_ids)
  if (is.null(features_are_rows)) {
    if (column_match && !row_match) {
      features_are_rows <- TRUE
    } else if (row_match && !column_match) {
      features_are_rows <- FALSE
    } else if (!row_match && !column_match) {
      stop("Count-table sample IDs must match metadata sample IDs exactly.",
           call. = FALSE)
    } else {
      stop("Could not infer count-table orientation uniquely; set `features_are_rows`.",
           call. = FALSE)
    }
  }
  if (!is.logical(features_are_rows) || length(features_are_rows) != 1L ||
      is.na(features_are_rows)) {
    stop("`features_are_rows` must be `TRUE`, `FALSE`, or `NULL`.", call. = FALSE)
  }
  if (!features_are_rows) counts <- t(counts)
  if (!setequal(colnames(counts), sample_ids)) {
    stop("Count-table sample IDs must match metadata sample IDs exactly.",
         call. = FALSE)
  }
  counts <- counts[, sample_ids, drop = FALSE]

  scores <- .ids_from_data_frame(
    decontam_result, feature_id_column, "decontam_result"
  )
  if (!score_column %in% names(scores) || !is.numeric(scores[[score_column]])) {
    stop("`decontam_result` must contain numeric column `", score_column, "`.",
         call. = FALSE)
  }
  if (any(scores[[score_column]] < 0 | scores[[score_column]] > 1, na.rm = TRUE)) {
    stop("Contaminant scores must be between 0 and 1 or `NA`.", call. = FALSE)
  }
  feature_ids <- rownames(counts)
  if (!setequal(rownames(scores), feature_ids)) {
    stop("Feature IDs in `decontam_result` and `count_table` must match exactly.",
         call. = FALSE)
  }
  scores <- scores[feature_ids, , drop = FALSE]

  if (!is.null(taxonomy)) {
    taxonomy <- .ids_from_data_frame(taxonomy, taxonomy_id_column, "taxonomy")
    if (!setequal(rownames(taxonomy), feature_ids)) {
      stop("Feature IDs in `taxonomy` and `count_table` must match exactly.",
           call. = FALSE)
    }
    taxonomy <- taxonomy[feature_ids, , drop = FALSE]
  }

  is_control <- metadata[[control_column]] %in% control_label
  if (!any(is_control) || all(is_control)) {
    stop("Both biological samples and negative controls are required.",
         call. = FALSE)
  }
  sample_group <- ifelse(is_control, "control", "biological")
  names(sample_group) <- sample_ids
  p <- scores[[score_column]]
  flags <- vapply(thresholds, function(th) !is.na(p) & p < th, logical(length(p)))
  if (is.null(dim(flags))) flags <- matrix(flags, ncol = 1L)
  rownames(flags) <- feature_ids
  colnames(flags) <- format(thresholds, trim = TRUE, scientific = FALSE)

  groups <- list(
    overall = rep(TRUE, length(sample_ids)),
    biological = !is_control,
    control = is_control
  )
  read_rows <- list()
  feature_rows <- list()
  sample_rows <- list()
  flag_rows <- list()

  for (i in seq_along(thresholds)) {
    retained <- !flags[, i]
    for (group_name in names(groups)) {
      sample_keep <- groups[[group_name]]
      before <- counts[, sample_keep, drop = FALSE]
      after <- counts[retained, sample_keep, drop = FALSE]
      reads_before <- sum(before)
      reads_after <- sum(after)
      features_before <- sum(rowSums(before) > 0)
      features_after <- sum(rowSums(after) > 0)
      read_rows[[length(read_rows) + 1L]] <- data.frame(
        threshold = thresholds[i], group = group_name,
        reads_before = reads_before, reads_retained = reads_after,
        reads_removed = reads_before - reads_after,
        reads_retained_pct = .safe_percent(reads_after, reads_before),
        stringsAsFactors = FALSE
      )
      feature_rows[[length(feature_rows) + 1L]] <- data.frame(
        threshold = thresholds[i], group = group_name,
        features_before = features_before, features_retained = features_after,
        features_removed = features_before - features_after,
        features_retained_pct = .safe_percent(features_after, features_before),
        stringsAsFactors = FALSE
      )
    }

    after_per_sample <- colSums(counts[retained, , drop = FALSE])
    before_per_sample <- colSums(counts)
    sample_rows[[i]] <- data.frame(
      threshold = thresholds[i], sample_id = sample_ids,
      group = sample_group, reads_before = before_per_sample,
      reads_retained = after_per_sample,
      reads_removed = before_per_sample - after_per_sample,
      reads_retained_pct = .safe_percent(after_per_sample, before_per_sample),
      zero_read_before = before_per_sample == 0,
      zero_read = after_per_sample == 0,
      stringsAsFactors = FALSE, row.names = NULL
    )
    flag_rows[[i]] <- data.frame(
      threshold = thresholds[i], feature_id = feature_ids,
      contaminant = flags[, i], stringsAsFactors = FALSE
    )
  }

  sample_totals <- colSums(counts)
  relative_counts <- sweep(counts, 2L, ifelse(sample_totals == 0, 1, sample_totals), "/")
  biological <- !is_control
  control <- is_control
  feature_summary <- data.frame(
    feature_id = feature_ids,
    score = p,
    reads_biological = rowSums(counts[, biological, drop = FALSE]),
    reads_control = rowSums(counts[, control, drop = FALSE]),
    mean_relative_abundance_biological = rowMeans(relative_counts[, biological, drop = FALSE]),
    mean_relative_abundance_control = rowMeans(relative_counts[, control, drop = FALSE]),
    prevalence_biological = rowSums(counts[, biological, drop = FALSE] > 0),
    prevalence_control = rowSums(counts[, control, drop = FALSE] > 0),
    prevalence_pct_biological = 100 * rowMeans(counts[, biological, drop = FALSE] > 0),
    prevalence_pct_control = 100 * rowMeans(counts[, control, drop = FALSE] > 0),
    stringsAsFactors = FALSE
  )
  first_flagged <- apply(flags, 1L, function(x) {
    hit <- which(x)
    if (length(hit)) thresholds[min(hit)] else NA_real_
  })
  feature_summary$first_flagged_threshold <- first_flagged
  if (!is.null(taxonomy)) {
    feature_summary <- cbind(feature_summary, taxonomy, row.names = NULL)
  }

  threshold_summary <- do.call(rbind, lapply(seq_along(thresholds), function(i) {
    ss <- sample_rows[[i]]
    data.frame(
      threshold = thresholds[i],
      contaminant_features = sum(flags[, i]),
      retained_features = sum(!flags[, i]),
      retained_samples = sum(!ss$zero_read),
      zero_read_samples = sum(ss$zero_read),
      biological_reads_retained_pct = read_rows[[3L * (i - 1L) + 2L]]$reads_retained_pct,
      control_reads_retained_pct = read_rows[[3L * (i - 1L) + 3L]]$reads_retained_pct,
      overall_reads_retained_pct = read_rows[[3L * (i - 1L) + 1L]]$reads_retained_pct,
      biological_features_retained_pct = feature_rows[[3L * (i - 1L) + 2L]]$features_retained_pct,
      control_features_retained_pct = feature_rows[[3L * (i - 1L) + 3L]]$features_retained_pct,
      overall_features_retained_pct = feature_rows[[3L * (i - 1L) + 1L]]$features_retained_pct
    )
  }))

  structure(list(
    thresholds = thresholds,
    threshold_summary = threshold_summary,
    read_retention = do.call(rbind, read_rows),
    feature_retention = do.call(rbind, feature_rows),
    sample_retention = do.call(rbind, sample_rows),
    feature_flags = do.call(rbind, flag_rows),
    feature_summary = feature_summary,
    scores = scores,
    counts = counts,
    metadata = metadata,
    taxonomy = taxonomy,
    is_control = stats::setNames(is_control, sample_ids),
    score_column = score_column,
    call = match.call()
  ), class = "decontam_sensitivity")
}

#' @export
print.decontam_sensitivity <- function(x, ...) {
  cat("decontam threshold sensitivity analysis\n")
  cat("  ", nrow(x$counts), " features; ", ncol(x$counts), " samples (",
      sum(x$is_control), " controls)\n", sep = "")
  print(x$threshold_summary, row.names = FALSE)
  invisible(x)
}
