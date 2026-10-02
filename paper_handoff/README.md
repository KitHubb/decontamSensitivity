# Skin microbiome decontamination study handoff

Open [`index.html`](index.html) for the complete PPT-oriented result guide and
the latest contamination-method benchmark reports.

## Final analysis decision

- Negative controls: Prep control + Air Swab
- Batch: DNA extraction Cycle aware
- Method: Frequency
- Threshold: 0.4
- Mock validation: Prep control only

## Contents

- `reports/`: self-contained HTML reports for basic analysis, decontamination, and mock validation
- `benchmarks/`: self-contained Test7 and Karstens benchmark reports, main figures, and compact result tables
- `figures/`: full-resolution PNG files selected for presentation
- `tables/`: compact CSV files supporting the reported values

## Latest benchmark reports

- `benchmarks/test7_actual_mock/04_2_Test7_actual_mock_benchmark.html`: Cell Mock and DNA Mock classification benchmark
- `benchmarks/test7_cross_method/06_Test7_combined_cross_method_validation.html`: decontam prevalence/frequency/combined, MicrobIEM, and SCRuB comparison
- `benchmarks/karstens2019/07_Karstens2019_three_level_benchmark.html`: Karstens et al. 2019 mock benchmark at lenient, default/middle, and aggressive settings

The preprocessing summary includes read depth and observed ASV counts for Total,
True samples, Controls, and Mocks. Both pooled observed ASVs and per-sample
Mean ± SD richness are provided.

The input phyloseq RDS is intentionally not tracked. To rerun the R Markdown files,
place `phy_F270R210_260921.rds` under `Phyloseq/` in the analysis working directory.
