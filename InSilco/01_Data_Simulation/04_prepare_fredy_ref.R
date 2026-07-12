# ==============================================================================
# Script: 04_prepare_fredy_ref.R
# Purpose: Convert TE annotations into FREDY-compatible BED4 and BED6 formats.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Teddy)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
})

data_dir <- "./data"
fredy_ref_dir <- "./results/FREDY_official_by_depth/ref"
dir.create(fredy_ref_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Load and process TE annotation
# ------------------------------------------------------------------------------
message("Loading TE annotation...")
te <- readRDS(file.path(data_dir, "mm10_TE.rds"))

if (exists("NCBI_check")) {
  te <- NCBI_check(te, ncbi_style = FALSE)
}

# ------------------------------------------------------------------------------
# 2. Build FREDY BED structures
# ------------------------------------------------------------------------------
message("Formatting BED files for FREDY...")
bed6 <- tibble(
  chrom = as.character(seqnames(te)),
  start = start(te) - 1L,
  end = end(te),
  name = as.character(mcols(te)$name),
  score = 0L,
  strand = as.character(strand(te))
) |>
  mutate(
    start = pmax(start, 0L),
    strand = if_else(strand %in% c("+", "-"), strand, "."),
    name = if_else(is.na(name) | name == "", as.character(mcols(te)$names), name),
    name = if_else(is.na(name) | name == "", paste0("TE_", row_number()), name)
  ) |>
  filter(!is.na(chrom), !is.na(start), !is.na(end), end > start) |>
  distinct() |>
  arrange(chrom, start, end, name)

bed4 <- bed6 |>
  select(chrom, start, end, name)

# ------------------------------------------------------------------------------
# 3. Export
# ------------------------------------------------------------------------------
write_tsv(bed4, file.path(fredy_ref_dir, "mm10_TE.FREDY.bed4"), col_names = FALSE)
write_tsv(bed6, file.path(fredy_ref_dir, "mm10_TE.FREDY.bed6"), col_names = FALSE)

message("FREDY references generated successfully in: ", fredy_ref_dir)