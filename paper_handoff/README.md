# Skin microbiome decontamination study handoff

Open [`index.html`](index.html) for the complete PPT-oriented result guide.

## Final analysis decision

- Negative controls: Prep control + Air Swab
- Batch: DNA extraction Cycle aware
- Method: Frequency
- Threshold: 0.4
- Mock validation: Prep control only

## Contents

- `reports/`: self-contained HTML reports for basic analysis, decontamination, and mock validation
- `figures/`: full-resolution PNG files selected for presentation
- `tables/`: compact CSV files supporting the reported values

The preprocessing summary includes read depth and observed ASV counts for Total,
True samples, Controls, and Mocks. Both pooled observed ASVs and per-sample
Mean ± SD richness are provided.

The input phyloseq RDS is intentionally not tracked. To rerun the R Markdown files,
place `phy_F270R210_260921.rds` under `Phyloseq/` in the analysis working directory.
