# ==============================================================================
# Script: 03_prepare_arriba_ref.R
# Purpose: Convert TE rds annotation into Arriba-compatible BED files.
# ==============================================================================
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
  library(readr)
})

data_dir <- "./data"
arriba_ref_dir <- "./results/arriba_merge/ref"
dir.create(arriba_ref_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load TE annotation
te <- readRDS(file.path(data_dir, "mm10_TE.rds"))

# 2. Build Arriba BED structure
bed8 <- tibble(
  chrom = as.character(seqnames(te)),
  start = start(te) - 1L,
  end = end(te),
  name = as.character(mcols(te)$name),
  score = 0L,
  strand = as.character(strand(te)),
  te_class = as.character(mcols(te)$class),
  te_family = as.character(mcols(te)$family)
) |>
  mutate(
    start = pmax(start, 0L),
    strand = if_else(strand %in% c("+", "-"), strand, "."),
    name = if_else(is.na(name) | name == "", as.character(mcols(te)$names), name),
    name = if_else(is.na(name) | name == "", paste0("TE_", row_number()), name),
    te_class = if_else(is.na(te_class) | te_class == "", "unknown", te_class),
    te_family = if_else(is.na(te_family) | te_family == "", "unknown", te_family)
  ) |>
  filter(!is.na(chrom), !is.na(start), !is.na(end), end > start) |>
  distinct() |>
  arrange(chrom, start, end, name)

# 3. Derive standard BED6 and BED4 for Arriba
bed6 <- bed8 |> select(chrom, start, end, name, score, strand)
bed4 <- bed8 |> select(chrom, start, end, name)

# 4. Export references
write_tsv(bed6, file.path(arriba_ref_dir, "mm10_TE.arriba.bed6"), col_names = FALSE)
write_tsv(bed4, file.path(arriba_ref_dir, "mm10_TE.arriba.bed4"), col_names = FALSE)

message("Arriba references generated successfully in: ", arriba_ref_dir)