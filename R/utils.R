.assert_scalar_character <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", name, "` must be one non-empty character value.", call. = FALSE)
  }
}

.ids_from_data_frame <- function(x, id_column, object_name) {
  if (!is.data.frame(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }

  if (!is.null(id_column)) {
    .assert_scalar_character(id_column, paste0(object_name, "_id_column"))
    if (!id_column %in% names(x)) {
      stop("Column `", id_column, "` is not present in `", object_name, "`.",
           call. = FALSE)
    }
    ids <- as.character(x[[id_column]])
    x[[id_column]] <- NULL
  } else {
    ids <- rownames(x)
  }

  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("`", object_name, "` must have unique, non-missing IDs.", call. = FALSE)
  }
  rownames(x) <- ids
  x
}

.validate_thresholds <- function(thresholds) {
  if (!is.numeric(thresholds) || !length(thresholds) || anyNA(thresholds) ||
      any(!is.finite(thresholds)) || any(thresholds < 0 | thresholds > 1)) {
    stop("`thresholds` must contain finite numeric values from 0 to 1.",
         call. = FALSE)
  }
  sort(unique(as.numeric(thresholds)))
}

.safe_percent <- function(numerator, denominator) {
  ifelse(denominator == 0, NA_real_, 100 * numerator / denominator)
}

.validate_group_colors <- function(group_colors) {
  required <- c("total", "biological", "control")
  if (!is.character(group_colors) || is.null(names(group_colors)) ||
      !all(required %in% names(group_colors)) ||
      anyNA(group_colors[required]) || any(!nzchar(group_colors[required]))) {
    stop("`group_colors` must be a named character vector containing `total`, `biological`, and `control`.",
         call. = FALSE)
  }
  group_colors[required]
}

.taxa_palette <- function(n) {
  if (n < 1L) return(character())
  paired <- RColorBrewer::brewer.pal(12L, "Paired")
  set1 <- RColorBrewer::brewer.pal(9L, "Set1")
  if (n <= 9L) return(paired[seq_len(n)])
  if (n <= 12L) return(grDevices::colorRampPalette(set1)(n))
  combined <- c(paired, set1)
  if (n <= 21L) return(combined[seq_len(n)])
  grDevices::colorRampPalette(combined)(n)
}

.qc_theme <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2),
      panel.grid.minor = ggplot2::element_line(linewidth = 0.1),
      text = ggplot2::element_text(size = 12),
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1)
    )
}

.threshold_index <- function(result, threshold) {
  if (!inherits(result, "decontam_sensitivity")) {
    stop("`result` must be returned by `run_threshold_sweep()`.", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
    stop("`threshold` must be one numeric value.", call. = FALSE)
  }
  index <- which(abs(result$thresholds - threshold) < sqrt(.Machine$double.eps))
  if (length(index) != 1L) {
    stop("Threshold ", threshold, " was not included in the sweep.", call. = FALSE)
  }
  index
}

.taxonomy_for_result <- function(result, taxonomy) {
  if (is.null(result$taxonomy)) {
    stop("No taxonomy was supplied to `run_threshold_sweep()`.", call. = FALSE)
  }
  .assert_scalar_character(taxonomy, "taxonomy")
  if (!taxonomy %in% names(result$taxonomy)) {
    stop("Taxonomy column `", taxonomy, "` is not available.", call. = FALSE)
  }
  values <- as.character(result$taxonomy[[taxonomy]])
  values[is.na(values) | !nzchar(values)] <- "Unclassified"
  values
}
