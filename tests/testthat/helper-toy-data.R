toy_sensitivity_data <- function(transpose = FALSE) {
  counts <- matrix(c(
    100, 80,  0,  0,
     10,  0, 20, 30,
      5,  5, 10,  0,
      0,  0, 15, 15,
     20, 20,  1,  1,
      2,  0,  0,  0
  ), nrow = 6, byrow = TRUE,
  dimnames = list(paste0("ASV", 1:6), c("B1", "B2", "C1", "C2")))
  metadata <- data.frame(
    type = c("sample", "sample", "control", "control"),
    row.names = colnames(counts)
  )
  scores <- data.frame(
    p = c(0.8, 0.05, 0.15, 0.25, 0.35, NA),
    row.names = rownames(counts)
  )
  taxonomy <- data.frame(
    Phylum = c(
      "Firmicutes", "Firmicutes", "Firmicutes",
      "Proteobacteria", "Firmicutes", NA
    ),
    Genus = c("A", "Staphylococcus", "Staphylococcus", "Ralstonia", "A", NA),
    row.names = rownames(counts)
  )
  list(
    counts = if (transpose) t(counts) else counts,
    metadata = metadata, scores = scores, taxonomy = taxonomy
  )
}

run_toy_sweep <- function(transpose = FALSE) {
  x <- toy_sensitivity_data(transpose)
  run_threshold_sweep(
    x$scores, x$counts, x$metadata, "type", "control",
    thresholds = c(0.1, 0.2, 0.4), taxonomy = x$taxonomy
  )
}
