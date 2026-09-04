#' Plot threshold sensitivity
#'
#' @inheritParams summarize_read_retention
#' @param metric One of `"both"`, `"reads"`, or `"features"`.
#' @return A `ggplot` object.
#' @export
plot_threshold_sensitivity <- function(result,
                                       metric = c("both", "reads", "features")) {
  if (!inherits(result, "decontam_sensitivity")) {
    stop("`result` must be returned by `run_threshold_sweep()`.", call. = FALSE)
  }
  metric <- match.arg(metric)
  reads <- result$read_retention[, c("threshold", "group", "reads_retained_pct")]
  names(reads)[3L] <- "retained_pct"
  reads$metric <- "Reads retained (%)"
  features <- result$feature_retention[, c("threshold", "group", "features_retained_pct")]
  names(features)[3L] <- "retained_pct"
  features$metric <- "Features retained (%)"
  plot_data <- switch(metric,
    reads = reads,
    features = features,
    both = rbind(reads, features)
  )
  plot_data$group <- factor(
    plot_data$group,
    levels = c("overall", "biological", "control"),
    labels = c("All samples", "Biological samples", "Negative controls")
  )
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = threshold, y = retained_pct, color = group, shape = group)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    ggplot2::scale_x_continuous(breaks = result$thresholds) +
    ggplot2::labs(x = "decontam threshold", y = NULL, color = NULL, shape = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
  if (metric == "both") {
    p <- p + ggplot2::facet_wrap(~metric, ncol = 1L)
  } else {
    p <- p + ggplot2::labs(y = unique(plot_data$metric))
  }
  p
}

#' Plot the distribution of decontam scores
#'
#' @inheritParams summarize_read_retention
#' @param bins Number of histogram bins.
#' @return A `ggplot` object with candidate thresholds marked.
#' @export
plot_decontam_scores <- function(result, bins = 50L) {
  if (!inherits(result, "decontam_sensitivity")) {
    stop("`result` must be returned by `run_threshold_sweep()`.", call. = FALSE)
  }
  ggplot2::ggplot(result$feature_summary, ggplot2::aes(x = score)) +
    ggplot2::geom_histogram(bins = bins, fill = "grey30", color = "white") +
    ggplot2::geom_vline(
      xintercept = result$thresholds, linetype = "dashed", color = "firebrick"
    ) +
    ggplot2::labs(x = "decontam score (p)", y = "Feature count") +
    ggplot2::theme_bw()
}
