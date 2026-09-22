#' Run a decontam threshold sweep directly from a phyloseq object
#'
#' Convenience wrapper for every contaminant-identification method currently
#' documented by `decontam`. Classification is obtained directly from
#' [decontam::isContaminant()] at every threshold instead of being inferred
#' from undocumented detailed-output columns. For `"either"` and `"both"`,
#' each sweep value is applied to both the frequency and prevalence tests.
#'
#' @param ps A `phyloseq` object.
#' @param control_column Sample-data column identifying negative controls.
#' @param control_label Value or values identifying negative controls.
#' @param thresholds Numeric thresholds to evaluate.
#' @param batch Optional sample-data column used as the `decontam` batch.
#' @param batch_combine Method used by `decontam` to combine batch scores.
#' @param normalize Passed to [decontam::isContaminant()].
#' @param method One of `"auto"`, `"frequency"`, `"prevalence"`,
#'   `"combined"`, `"minimum"`, `"either"`, or `"both"`.
#' @param concentration_column Sample-data column containing positive numeric
#'   DNA concentrations. Required by every method with a frequency component.
#' @return A `decontam_sensitivity` object.
#' @export
run_decontam_threshold_sweep <- function(ps,
                                         control_column,
                                         control_label,
                                         thresholds = seq(0.1, 0.5, 0.1),
                                         batch = NULL,
                                         batch_combine = c("minimum", "product", "fisher"),
                                         normalize = TRUE,
                                         method = c(
                                           "prevalence", "frequency", "combined",
                                           "minimum", "either", "both", "auto"
                                         ),
                                         concentration_column = NULL) {
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package `phyloseq` is required for this function.", call. = FALSE)
  }
  if (!requireNamespace("decontam", quietly = TRUE)) {
    stop("Package `decontam` is required for this function.", call. = FALSE)
  }
  thresholds <- .validate_thresholds(thresholds)
  method <- match.arg(method)
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

  frequency_methods <- c("frequency", "combined", "minimum", "either", "both")
  needs_concentration <- method %in% frequency_methods ||
    identical(method, "auto") && !is.null(concentration_column)
  concentration <- NULL
  if (needs_concentration) {
    if (is.null(concentration_column)) {
      stop(
        "`concentration_column` is required for method `", method, "`.",
        call. = FALSE
      )
    }
    .assert_scalar_character(concentration_column, "concentration_column")
    if (!concentration_column %in% names(metadata)) {
      stop(
        "Column `", concentration_column, "` is not present in sample data.",
        call. = FALSE
      )
    }
    concentration <- metadata[[concentration_column]]
    if (!is.numeric(concentration) || anyNA(concentration) ||
        any(!is.finite(concentration)) || any(concentration <= 0)) {
      stop(
        "`concentration_column` must contain finite, positive numeric values.",
        call. = FALSE
      )
    }
  }

  effective_method <- if (identical(method, "auto")) {
    if (is.null(concentration)) "prevalence" else "combined"
  } else method

  call_decontam <- function(threshold, detailed) {
    threshold_arg <- if (effective_method %in% c("either", "both")) {
      c(threshold, threshold)
    } else threshold
    args <- list(
      seqtab = ps,
      method = effective_method,
      batch = batch_values,
      batch.combine = batch_combine,
      threshold = threshold_arg,
      normalize = normalize,
      detailed = detailed
    )
    if (identical(effective_method, "frequency")) {
      args$conc <- concentration
    } else if (identical(effective_method, "prevalence")) {
      args$neg <- is_control
    } else {
      args$conc <- concentration
      args$neg <- is_control
    }
    do.call(decontam::isContaminant, args)
  }

  scores <- call_decontam(thresholds[1L], detailed = TRUE)
  if (!is.data.frame(scores)) {
    stop(
      "`decontam::isContaminant(..., detailed = TRUE)` did not return a data frame; check the installed decontam version.",
      call. = FALSE
    )
  }

  feature_ids <- phyloseq::taxa_names(ps)
  if (is.null(rownames(scores)) || !setequal(rownames(scores), feature_ids)) {
    stop(
      "Detailed decontam output did not contain the expected feature IDs.",
      call. = FALSE
    )
  }
  scores <- scores[feature_ids, , drop = FALSE]

  flags <- vapply(thresholds, function(threshold) {
    out <- call_decontam(threshold, detailed = FALSE)
    if (!is.logical(out) || length(out) != length(feature_ids) || anyNA(out)) {
      stop(
        "`decontam::isContaminant(..., detailed = FALSE)` returned an incompatible classification vector.",
        call. = FALSE
      )
    }
    if (!is.null(names(out))) {
      if (!setequal(names(out), feature_ids)) {
        stop(
          "decontam classifications did not contain the expected feature IDs.",
          call. = FALSE
        )
      }
      out <- out[feature_ids]
    }
    unname(out)
  }, logical(length(feature_ids)))
  if (is.null(dim(flags))) flags <- matrix(flags, ncol = 1L)
  rownames(flags) <- feature_ids

  numeric_column <- function(name) {
    if (name %in% names(scores) && is.numeric(scores[[name]])) {
      scores[[name]]
    } else {
      rep(NA_real_, nrow(scores))
    }
  }
  overall_p <- numeric_column("p")
  frequency_p <- numeric_column("p.freq")
  prevalence_p <- numeric_column("p.prev")
  sensitivity_score <- switch(
    effective_method,
    frequency = if (all(is.na(overall_p))) frequency_p else overall_p,
    prevalence = if (all(is.na(overall_p))) prevalence_p else overall_p,
    combined = overall_p,
    minimum = overall_p,
    either = {
      value <- pmin(frequency_p, prevalence_p, na.rm = TRUE)
      value[is.infinite(value)] <- NA_real_
      value
    },
    both = pmax(frequency_p, prevalence_p)
  )
  score_definition <- switch(
    effective_method,
    either = "min(p.freq, p.prev)",
    both = "max(p.freq, p.prev)",
    "decontam p"
  )
  score_matches_flags <- vapply(seq_along(thresholds), function(i) {
    expected <- !is.na(sensitivity_score) & sensitivity_score < thresholds[i]
    identical(unname(expected), unname(flags[, i]))
  }, logical(1))
  if (!all(score_matches_flags)) {
    sensitivity_score[] <- NA_real_
    score_definition <- "unavailable; classifications use direct decontam calls"
    warning(
      "Detailed decontam scores were not classification-equivalent to the direct results; the score plot has been disabled, but threshold classifications remain valid.",
      call. = FALSE
    )
  }
  scores$.sensitivity_score <- sensitivity_score

  counts <- .phyloseq_count_matrix(ps)
  taxonomy <- NULL
  if (!is.null(phyloseq::tax_table(ps, errorIfNULL = FALSE))) {
    taxonomy <- as.data.frame(phyloseq::tax_table(ps), stringsAsFactors = FALSE)
  }
  result <- run_threshold_sweep(
    decontam_result = scores, count_table = counts, metadata = metadata,
    control_column = control_column, control_label = control_label,
    thresholds = thresholds, taxonomy = taxonomy,
    score_column = ".sensitivity_score", contaminant_flags = flags,
    features_are_rows = TRUE
  )
  result$method <- effective_method
  result$requested_method <- method
  result$concentration_column <- concentration_column
  result$score_definition <- score_definition
  result$score_available <- all(score_matches_flags)
  result$score_label <- if (!result$score_available) {
    "Scalar score unavailable; direct decontam classifications used"
  } else if (effective_method %in% c("either", "both")) {
    paste0("Threshold-equivalent score: ", score_definition)
  } else {
    paste0("decontam ", effective_method, " score (p)")
  }
  result$decontam_version <- as.character(utils::packageVersion("decontam"))
  result$call <- match.call()
  result
}
