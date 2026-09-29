# Test9 original-paper figure reanalysis

This handoff contains the self-contained HTML report, the source R Markdown,
major PNG figures, and compact result tables from script 13.

## Analysis definition

- Input: Test9 DADA2/GG2 phyloseq produced by script 09.
- Contaminant removal: pooled `decontam` prevalence method, threshold 0.5,
  matching script 12 (skin swabs versus negative extraction controls).
- Clinical population: the 2,319 original Supplementary Table S1 biological
  samples, representing 418 individuals. Technical SRA runs mapped to the same
  original sample were summed.
- Alpha/beta diversity: rarefied to 5,000 reads with seed 42. All 2,319 samples
  remained eligible.
- Taxonomic composition: unrarefied decontaminated counts converted to
  within-sample relative abundance. Colors come from the existing
  `functions.R::taxa_plot()` phylum palette.

## Main draft result

The original lesion-gradient pattern was recovered: corresponding control skin
had the highest Shannon diversity, contralateral patient skin was intermediate,
and perilesional patient skin had the lowest diversity. The three planned
Shannon contrasts remained significant after BH correction.

Bray-Curtis and Jaccard PERMANOVA results were significant for all displayed
clinical comparisons, although effect sizes were small. Upper-back Bray-Curtis
also showed significant dispersion heterogeneity, so that comparison cannot be
interpreted solely as centroid separation.

## Limits to review

- The PCoA plots are unconstrained descriptive ordinations; they are not the
  original paper's partial constrained ordinations.
- The current GG2 taxonomy contains no ASV classified exactly as
  `Staphylococcus aureus`. The BPDAI table matched 210 of 228 patients, but the
  BPDAI-S. aureus correlation was not estimated rather than substituting another
  Staphylococcus taxon.
- The Supplementary Table S1 source includes `Extraction13` and the apparent
  typo `Extraxtion13` as separate labels. This raw-label issue is visible in the
  batch table and should be resolved before treating extraction-batch estimates
  as manuscript results.
- Several batches contain only one status group or one center. Their batch
  effects cannot be separated from the corresponding clinical/site structure.

Open `13_Test9_original_paper_replication.html` for the complete report.
