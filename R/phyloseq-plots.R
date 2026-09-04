#' Plot biological-versus-control prevalence enrichment
#'
#' Replaces the complete `kA`, `odds.sample`, `odds.control`, `p1`, `p2`, and
#' `ggarrange()` block with one function call.
#'
#' @inheritParams summarize_control_prevalence
#' @inheritParams calculate_prevalence_enrichment
#' @param log10_x Logical; use a log10 scale for the prevalence ratio.
#' @return A ggplot object. The underlying feature table is available through
#'   `summarize_control_prevalence()`.
#' @export
plot_prevalence_enrichment <- function(ps,
                                       control_column,
                                       control_label,
                                       prevalence_unit = c("count", "proportion"),
                                       pseudocount = 0,
                                       log10_x = FALSE) {
  prevalence_unit <- match.arg(prevalence_unit)
  summary <- calculate_prevalence_enrichment(
    ps = ps,
    control_column = control_column,
    control_label = control_label,
    prevalence_unit = prevalence_unit,
    pseudocount = pseudocount
  )
  plot_data <- rbind(
    data.frame(
      enrichment = summary$odds.sample,
      abundance_total = summary$summarized,
      direction = "a) Enrichment in Samples"
    ),
    data.frame(
      enrichment = summary$odds.control,
      abundance_total = summary$summarized,
      direction = "b) Enrichment in Controls"
    )
  )
  plot_data <- plot_data[
    is.finite(plot_data$enrichment) & plot_data$enrichment > 0 &
      plot_data$abundance_total > 0,
    ,
    drop = FALSE
  ]
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = enrichment, y = abundance_total)
  ) +
    ggplot2::geom_point(alpha = 0.65) +
    ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~direction, ncol = 2, scales = "free_x") +
    ggplot2::labs(
      x = paste0("Prevalence enrichment ratio (", prevalence_unit, ")"),
      y = "taxa_sum"
    ) +
    ggplot2::theme_bw()
  if (isTRUE(log10_x)) p <- p + ggplot2::scale_x_log10()
  p
}

#' Plot sample-level read retention as boxplots
#'
#' @inheritParams summarize_read_retention
#' @param show_points Logical; overlay individual samples.
#' @return A ggplot object.
#' @export
plot_sample_read_retention <- function(result, show_points = TRUE) {
  if (!inherits(result, "decontam_sensitivity")) {
    stop("`result` must be returned by `run_threshold_sweep()`.", call. = FALSE)
  }
  plot_data <- result$sample_retention
  plot_data$threshold_label <- factor(
    format(plot_data$threshold, trim = TRUE, scientific = FALSE),
    levels = format(result$thresholds, trim = TRUE, scientific = FALSE)
  )
  plot_data$group <- factor(
    plot_data$group,
    levels = c("biological", "control"),
    labels = c("Biological samples", "Negative controls")
  )
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = threshold_label, y = reads_retained_pct, fill = group)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    ggplot2::labs(
      x = "decontam threshold",
      y = "Reads retained per sample (%)",
      fill = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
  if (isTRUE(show_points)) {
    p <- p + ggplot2::geom_point(
      position = ggplot2::position_jitterdodge(
        jitter.width = 0.12, dodge.width = 0.75
      ),
      alpha = 0.45,
      size = 1
    )
  }
  p
}

#' Plot raw reads assigned to flagged taxa
#'
#' @inheritParams summarize_flagged_taxa
#' @param top_n Maximum number of taxa to display.
#' @param measure Either `"mean_per_sample"` (recommended when groups contain
#'   different sample counts) or `"total"`.
#' @return A ggplot object comparing biological samples and controls.
#' @export
plot_flagged_taxa_reads <- function(result,
                                    threshold,
                                    taxonomy = "Genus",
                                    top_n = 20L,
                                    measure = c("mean_per_sample", "total")) {
  measure <- match.arg(measure)
  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 1) {
    stop("`top_n` must be one positive number.", call. = FALSE)
  }
  taxa <- summarize_flagged_taxa(result, threshold, taxonomy)
  if (!nrow(taxa)) {
    stop("No taxa are flagged at the requested threshold.", call. = FALSE)
  }
  taxa$total_flagged_reads <- taxa$reads_biological + taxa$reads_control
  taxa <- utils::head(
    taxa[order(-taxa$total_flagged_reads), , drop = FALSE], top_n
  )
  n_biological <- sum(!result$is_control)
  n_control <- sum(result$is_control)
  if (measure == "mean_per_sample") {
    biological_values <- taxa$reads_biological / n_biological
    control_values <- taxa$reads_control / n_control
    y_label <- "Mean raw reads per sample"
  } else {
    biological_values <- taxa$reads_biological
    control_values <- taxa$reads_control
    y_label <- "Total raw read count"
  }
  plot_data <- rbind(
    data.frame(
      taxon = taxa$taxon, group = "Biological samples",
      reads = biological_values
    ),
    data.frame(
      taxon = taxa$taxon, group = "Negative controls",
      reads = control_values
    )
  )
  plot_data$taxon <- factor(plot_data$taxon, levels = rev(taxa$taxon))
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = taxon, y = reads, fill = group)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = taxonomy,
      y = y_label,
      fill = NULL,
      title = paste("Reads assigned to flagged", taxonomy, "at threshold", threshold)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}

#' Plot taxon raw-read composition before and after decontamination
#'
#' @param ps A phyloseq object used to produce `result`.
#' @inheritParams summarize_flagged_taxa
#' @inheritParams split_phyloseq_groups
#' @param top_n Number of most abundant taxa to display individually.
#' @param other_label Label used for all remaining taxa.
#' @return A faceted sample-level raw-read stacked-bar ggplot.
#' @export
plot_taxa_reads_before_after <- function(ps,
                                         result,
                                         threshold,
                                         control_column,
                                         control_label,
                                         taxonomy = "Genus",
                                         top_n = 15L,
                                         other_label = "Other") {
  group_info <- .phyloseq_group_info(ps, control_column, control_label)
  .assert_scalar_character(taxonomy, "taxonomy")
  if (!taxonomy %in% phyloseq::rank_names(ps)) {
    stop("Taxonomic rank `", taxonomy, "` is not available in `ps`.",
         call. = FALSE)
  }
  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 1) {
    stop("`top_n` must be one positive number.", call. = FALSE)
  }
  filtered <- filter_phyloseq_at_threshold(ps, result, threshold)

  melt_taxa <- function(object, state) {
    glommed <- phyloseq::tax_glom(object, taxrank = taxonomy, NArm = FALSE)
    out <- phyloseq::psmelt(glommed)
    out$taxon <- as.character(out[[taxonomy]])
    out$taxon[is.na(out$taxon) | !nzchar(out$taxon)] <- "Unclassified"
    out$state <- state
    out
  }
  before_label <- "Before decontam"
  after_label <- paste("After threshold", threshold)
  before <- melt_taxa(ps, before_label)
  after <- melt_taxa(filtered, after_label)
  totals <- stats::aggregate(Abundance ~ taxon, data = before, FUN = sum)
  totals <- totals[order(-totals$Abundance), , drop = FALSE]
  top_taxa <- utils::head(totals$taxon, top_n)
  plot_data <- rbind(before, after)
  plot_data$taxon_plot <- ifelse(
    plot_data$taxon %in% top_taxa, plot_data$taxon, other_label
  )
  plot_data$state <- factor(
    plot_data$state, levels = c(before_label, after_label)
  )
  sample_control <- group_info$is_control[as.character(plot_data$Sample)]
  plot_data$sample_group <- factor(
    ifelse(sample_control, "control", "biological"),
    levels = c("biological", "control"),
    labels = c("Biological samples", "Negative controls")
  )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = Sample, y = Abundance, fill = taxon_plot)
  ) +
    ggplot2::geom_col(width = 1) +
    ggplot2::facet_grid(
      state ~ sample_group, scales = "free_x", space = "free_x"
    ) +
    ggplot2::labs(x = NULL, y = "Raw read count", fill = taxonomy) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}
