# decontamSensitivity

`decontamSensitivity` provides QC summaries and plots for comparing
`decontam` thresholds in microbiome data.

## Workspace layout

- `_paper/Review/`: literature searches, screening evidence, reference checks, and review tools (local only).
- `_paper/Docs/`: workspace documentation (local only).
- `_paper/Article/`: manuscript drafts, editorial notes, and Word preparation tools (local only).
- `Analysis_R/`: analysis examples, data preparation scripts, and analysis results.

Open `decontamSensitivity.Rproj` in this root directory. The package directories
`R/`, `man/`, `tests/`, `vignettes/`, and `data/` remain here for R package compatibility.
Run analysis scripts with the project root as the working directory.

## Why decontamSensitivity?

Choosing a `decontam` threshold can be difficult in low-biomass studies because
contaminant removal may also remove biological signal.

`decontamSensitivity` helps evaluate this trade-off across multiple thresholds by asking:

- How many reads and features are retained in biological samples?
- How strongly are reads and features reduced in negative controls?
- Which taxa are affected as the threshold changes?

The package does not choose an “optimal” threshold. It provides evidence for a
choice based on the controls, expected biology, and study design.

## Installation

Install the latest version from GitHub:

```r
install.packages("pak") # Skip this if pak is already installed
pak::pak("KitHubb/decontamSensitivity")

library(decontamSensitivity)
```

## Quick Start

Run the full workflow with one function. The default method is prevalence:

`ps` must be a `phyloseq` object containing an OTU table and sample metadata.
The metadata column supplied to `control_column` must identify the negative
controls. OTU tables work in either orientation.

```r
qc <- run_decontam_qc(
  ps = ps,
  control_column = "sample_type",
  control_label = "control",
  thresholds = seq(0.1, 0.9, by = 0.1),
  taxonomy = "Genus",
  group_colors = c(
    total = "black",
    biological = "tomato",
    control = "steelblue"
  )
)
```

```r
# Summary tables
qc$tables$threshold_summary
qc$tables$sample_retention_summary
qc$tables$flagged_taxa

# Main plots
qc$plots$threshold_sensitivity
qc$plots$sample_retention
qc$plots$prevalence_enrichment

# Results at threshold 0.5
qc$plots$flagged_taxa_reads_by_threshold[["0.5"]]
qc$plots$taxa_reads_before_after_by_threshold[["0.5"]]
ps_filtered <- qc$filtered_phyloseq_by_threshold[["0.5"]]
```

Taxon colors are chosen automatically from `RColorBrewer`. Use `taxa_colors`
to set your own.

## Other Functions

Use each step separately when you want more control.

```r
result <- run_decontam_threshold_sweep(
  ps = ps,
  control_column = "sample_type",
  control_label = "control",
  thresholds = seq(0.1, 0.9, by = 0.1)
)

summarize_read_retention(result)
summarize_feature_retention(result)
summarize_sample_read_retention(result)
summarize_flagged_taxa(result, threshold = 0.5, taxonomy = "Genus")

plot_decontam_scores(result)
plot_threshold_sensitivity(result)
plot_sample_read_retention(result)
plot_flagged_taxa_reads(result, threshold = 0.5, taxonomy = "Genus")
```

Filter the original object or split biological samples and controls:

```r
ps_filtered <- filter_phyloseq_at_threshold(ps, result, threshold = 0.5)

groups <- split_phyloseq_groups(
  ps,
  control_column = "sample_type",
  control_label = "control"
)
```

Check whether taxa are more common in biological samples or controls:

```r
enrichment <- calculate_prevalence_enrichment(
  ps,
  control_column = "sample_type",
  control_label = "control"
)

plot_prevalence_enrichment(
  ps,
  control_column = "sample_type",
  control_label = "control"
)
```

The columns `odds.sample` and `odds.control` are prevalence ratios, not odds
ratios.

### Supported decontam methods

`run_decontam_threshold_sweep()` and `run_decontam_qc()` support every method
documented by `decontam`: `"auto"`, `"frequency"`, `"prevalence"`,
`"combined"`, `"minimum"`, `"either"`, and `"both"`. The official method name
is `"minimum"`, not `"minimal"`.

Methods with a frequency component require a positive numeric concentration
column in `sample_data(ps)`:

```r
frequency_result <- run_decontam_threshold_sweep(
  ps,
  control_column = "sample_type",
  control_label = "control",
  thresholds = seq(0.1, 0.9, by = 0.1),
  method = "frequency",
  concentration_column = "DNA_concentration"
)

both_result <- run_decontam_threshold_sweep(
  ps,
  control_column = "sample_type",
  control_label = "control",
  thresholds = seq(0.1, 0.9, by = 0.1),
  method = "both",
  concentration_column = "DNA_concentration"
)
```

For `"either"` and `"both"`, each value in `thresholds` is applied to both the
frequency and prevalence tests. Classification flags come from direct
`decontam::isContaminant(..., detailed = FALSE)` calls at every threshold;
detailed-output score columns are used only for diagnostic plots. See the
[all-methods example](vignettes/hv-threshold-sensitivity.Rmd#all-decontam-methods)
for more examples.

### Test-dataset method comparison

The repository includes a reproducible comparison of the frequency,
prevalence, and combined methods on the small test dataset. It evaluates
thresholds `0.01`, `0.05`, and `0.1` through `0.9` in increments of `0.1`,
then exports threshold summaries, feature-level flags, pairwise method
agreement, and a comparison plot:

```r
source("Analysis_R/examples/test_dataset_method_threshold_comparison.R")
```

Outputs are written to the analysis workspace at
`../decontamSensitivity_paper/R_analysis/Results/test_method_comparison/`.
Set the `DECONTAM_ANALYSIS_RESULTS` environment variable only when that
workspace is stored elsewhere. This toy dataset is suitable for checking code
behavior, but it is too small to justify an optimal threshold for a real study.

## Published use case

These plots come from a published healthy-volunteer skin microbiome study.

### QC overview

![HV decontam QC overview](man/figures/hv_qc_overview.png)

The package plots can be combined into one figure with `patchwork`:

```r
library(patchwork)
library(ggplot2)

p_threshold <- qc$plots$threshold_sensitivity
p_flagged <- qc$plots$flagged_taxa_reads_by_threshold[["0.1"]]
p_taxa <- qc$plots$taxa_reads_before_after_by_threshold[["0.1"]] +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

qc_figure <-
  (p_threshold | p_flagged) /
  p_taxa +
  plot_layout(heights = c(1, 1.15), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

qc_figure
```

### How the threshold was chosen

In this study, only three negative controls were available, so the threshold
was interpreted cautiously rather than selected from a single metric.
Thresholds of 0.1–0.3 produced similar overall results, whereas filtering at
0.4 removed *Cutibacterium*, a common member of the skin microbiome. We
therefore considered 0.1–0.3 to be a reasonable range. A threshold of 0.1 was
selected because it gave results comparable to 0.3 while taking the more
conservative approach of preserving as much biological signal as possible.

### Threshold sensitivity

![Threshold-specific read and feature retention](man/figures/hv_threshold_sensitivity.png)

### Sample read retention

![Sample-level read retention](man/figures/hv_sample_retention.png)

### Flagged taxa

![Flagged genus read counts](man/figures/hv_flagged_taxa_threshold_0.1.png)

![Flagged genus read counts across thresholds](man/figures/hv_flagged_taxa_by_threshold.gif)

### Taxa composition

![Taxa composition before and after filtering](man/figures/hv_taxa_reads_before_after_by_threshold.gif)

More details: [published HV case study](https://www.nature.com/articles/s41598-026-62903-7).

## Citation

Please cite `decontam`, which provides the contaminant-identification method:

> Davis NM, Proctor DM, Holmes SP, Relman DA, Callahan BJ. Simple statistical
> identification and removal of contaminant sequences in marker-gene and
> metagenomics data. *Microbiome*. 2018;6:226.
> [https://doi.org/10.1186/s40168-018-0605-2](https://doi.org/10.1186/s40168-018-0605-2)

```r
citation("decontam")
```
