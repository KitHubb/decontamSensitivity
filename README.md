# decontamSensitivity

`decontamSensitivity` is a companion utility for evaluating threshold
sensitivity of contaminant classifications generated with the R package
[`decontam`](https://github.com/benjjneb/decontam). It does not reimplement
the contaminant-identification algorithms in `decontam`, and it is not
affiliated with or endorsed by the `decontam` developers.

The package answers a practical QC question: as the classification threshold
changes, how much biological signal is preserved and how much negative-control
signal is depleted? It reports:

- read retention in biological samples, negative controls, and all samples;
- feature retention using group-specific observed-feature denominators;
- retained and zero-read samples at each threshold;
- features and taxa newly flagged between thresholds; and
- the biological-sample abundance of newly flagged taxa.

It is a sensitivity diagnostic, not an automatic threshold optimizer.

## Installation

Install the latest version directly from GitHub with
[`pak`](https://pak.r-lib.org/):

```r
install.packages("pak") # Run once if pak is not installed
pak::pak("KitHubb/decontamSensitivity")

library(decontamSensitivity)
packageVersion("decontamSensitivity") # 0.1.1 or later
```

`pak` installs the required `decontam`, `phyloseq`, and `ggplot2`
dependencies declared by the package. Restart the R session before reinstalling
an already-loaded version, then verify that the public interface is available:

```r
packageVersion("decontamSensitivity")
"run_decontam_qc" %in% getNamespaceExports("decontamSensitivity")
```

## Citation

`decontamSensitivity` is a diagnostic companion to `decontam` and does not
replace its contaminant-identification method. Analyses using this package
should therefore cite the original `decontam` publication:

> Davis NM, Proctor DM, Holmes SP, Relman DA, Callahan BJ. Simple statistical
> identification and removal of contaminant sequences in marker-gene and
> metagenomics data. *Microbiome*. 2018;6:226.
> [https://doi.org/10.1186/s40168-018-0605-2](https://doi.org/10.1186/s40168-018-0605-2)

```bibtex
@article{Davis2018decontam,
  author  = {Davis, Nicole M. and Proctor, Diana M. and Holmes, Susan P. and
             Relman, David A. and Callahan, Benjamin J.},
  title   = {Simple statistical identification and removal of contaminant
             sequences in marker-gene and metagenomics data},
  journal = {Microbiome},
  year    = {2018},
  volume  = {6},
  pages   = {226},
  doi     = {10.1186/s40168-018-0605-2}
}
```

The citation bundled with the installed upstream package can also be checked
from R:

```r
citation("decontam")
```

## One-call phyloseq workflow

For an existing phyloseq object, one call runs the complete QC workflow. The
selected threshold is always supplied by the analyst; the package does not
optimize it automatically.

```r
qc <- run_decontam_qc(
  ps,
  control_column = "T_C",
  control_label = "Ctrl",
  thresholds = seq(0.1, 0.5, 0.1),
  selected_threshold = 0.1,
  batch = "SequencingRun",
  plot_all_thresholds = TRUE
)

qc$tables$threshold_summary
qc$tables$library_size_summary
qc$tables$prevalence_comparison
qc$tables$prevalence_enrichment
qc$tables$flagged_taxa
qc$tables$newly_flagged_taxa

qc$plots$sample_retention
qc$plots$prevalence_enrichment
qc$plots$flagged_taxa_reads
qc$plots$taxa_reads_before_after
qc$plots$flagged_taxa_reads_by_threshold[["0.5"]]
qc$plots$taxa_reads_before_after_by_threshold[["0.5"]]

ps_filtered <- qc$filtered_phyloseq
```

The plot and summary functions are also exported individually when only one
component is needed, including `plot_sample_read_retention()`,
`plot_flagged_taxa_reads()`, `plot_taxa_reads_before_after()`,
`calculate_prevalence_enrichment()`, `plot_prevalence_enrichment()`,
`summarize_control_prevalence()`, and `split_phyloseq_groups()`.

The original sample/control enrichment calculation is available directly:

```r
kA <- calculate_prevalence_enrichment(
  ps,
  control_column = "T_C",
  control_label = "Ctrl",
  prevalence_unit = "count",
  pseudocount = 0
)

kA$odds.sample
kA$odds.control

plot_prevalence_enrichment(
  ps,
  control_column = "T_C",
  control_label = "Ctrl"
)
```

These legacy `odds.*` columns are prevalence ratios, matching the original
analysis code; they are not statistical odds ratios.

## Score-table workflow

Already-calculated `decontam` scores can also be evaluated directly. No study
or example dataset is distributed with this repository.

```r
library(decontamSensitivity)

result <- run_threshold_sweep(
  decontam_result = decontam_scores,
  count_table = feature_counts,
  metadata = sample_metadata,
  control_column = "sample_class",
  control_label = "control",
  thresholds = seq(0.1, 0.5, 0.1),
  taxonomy = feature_taxonomy
)

result$threshold_summary
plot_threshold_sensitivity(result)
plot_decontam_scores(result)

summarize_flagged_taxa(result, threshold = 0.1, taxonomy = "Genus")
compare_newly_flagged_taxa(
  result,
  threshold_low = 0.1,
  threshold_high = 0.4,
  taxonomy = "Genus"
)
```

See [`examples/published_skin_microbiome_example.R`](examples/published_skin_microbiome_example.R)
for the schema used by the published workflow.

The DT Swab dataset can be evaluated directly with the package using
[`examples/DT_swab_threshold_sensitivity.Rmd`](examples/DT_swab_threshold_sensitivity.Rmd).
This report restores the five negative-control labels from the study metadata,
runs thresholds 0.1--0.9, summarizes newly flagged taxa, and validates the
threshold 0.5 result against the previously saved v4 phyloseq object.

## Interpretation

A threshold can be considered conservative when negative-control reads are
substantially depleted while biological reads and biologically plausible taxa
remain stable. Inspect newly flagged taxa whenever a higher threshold causes a
large loss. Beta-diversity separation between biological samples and controls
is deliberately not used as an optimization criterion because the prevalence
model already uses control status.

## Published use case

The initial workflow was used in the threshold sensitivity analysis reported
in [Scientific Reports (2026)](https://www.nature.com/articles/s41598-026-62903-7#Sec17).
The original `decontam` method listed in the [Citation](#citation) section
should be cited alongside this workflow.
