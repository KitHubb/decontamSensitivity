#' Run a prevalence sweep directly from a phyloseq object
#'
#' Convenience wrapper that calculates one set of `decontam` prevalence
#' scores, then passes those scores and the aligned data to
#' [run_threshold_sweep()]. Thresholds are applied downstream, so the
#' statistical model is not refit for every threshold.
#'
#' @param ps A `phyloseq` object.
#' @param control_column Sample-data column identifying negative controls.
#' @param control_label Value or values identifying negative controls.
#' @param thresholds Numeric thresholds to evaluate.
#' @param batch Optional sample-data column used as the `decontam` batch.
#' @param batch_combine Method used by `decontam` to combine batch scores.
#' @param normalize Passed to [decontam::isContaminant()].
#' @return A `decontam_sensitivity` object.
#' @export
run_decontam_threshold_sweep <- function(ps,
                                         control_column,
                                         control_label,
                                         thresholds = seq(0.1, 0.5, 0.1),
                                         batch = NULL,
                                         batch_combine = c("minimum", "product", "fisher"),
                                         normalize = TRUE) {
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package `phyloseq` is required for this function.", call. = FALSE)
  }
  if (!requireNamespace("decontam", quietly = TRUE)) {
    stop("Package `decontam` is required for this function.", call. = FALSE)
  }
  thresholds <- .validate_thresholds(thresholds)
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
  is_control <- metadata[[control_column]] %in% control_label
  if (!any(is_control) || all(is_control)) {
    stop("Both biological samples and negative controls are required.",
         call. = FALSE)
  }
  batch_combine <- match.arg(batch_combine)
  batch_values <- NULL
  if (!is.null(batch)) {
    .assert_scalar_character(batch, "batch")
    if (!batch %in% names(metadata)) {
      stop("Column `", batch, "` is not present in sample data.", call. = FALSE)
    }
    batch_values <- metadata[[batch]]
  }

  scores <- decontam::isContaminant(
    ps, method = "prevalence", neg = is_control, batch = batch_values,
    batch.combine = batch_combine, threshold = thresholds[1L],
    normalize = normalize, detailed = TRUE
  )
  counts <- .phyloseq_count_matrix(ps)
  taxonomy <- NULL
  if (!is.null(phyloseq::tax_table(ps, errorIfNULL = FALSE))) {
    taxonomy <- as.data.frame(phyloseq::tax_table(ps), stringsAsFactors = FALSE)
  }
  run_threshold_sweep(
    decontam_result = scores, count_table = counts, metadata = metadata,
    control_column = control_column, control_label = control_label,
    thresholds = thresholds, taxonomy = taxonomy, features_are_rows = TRUE
  )
}
