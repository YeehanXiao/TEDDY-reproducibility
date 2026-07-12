# ==============================================================================
# Script: 05_prepare_lions_ref.R
# Purpose: Format TE annotation into a 7-column BED without 'chr' prefix
#          specifically required by LIONS.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
  library(readr)
})

data_dir <- "./data"
lions_ref_dir <- "./results/LIONS_official_by_depth/ref"
dir.create(lions_ref_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Load TE Annotation
# ------------------------------------------------------------------------------
te_path <- file.path(data_dir, "mm10_TE.rds")
if (!file.exists(te_path)) stop("Missing TE annotation: ", te_path)

te <- readRDS(te_path)

# ------------------------------------------------------------------------------
# 2. Format for LIONS (7 columns, no 'chr')
# ------------------------------------------------------------------------------
mc <- mcols(te)

name <- if ("name" %in% colnames(mc)) as.character(mc$name) else as.character(mc$names)
cls  <- if ("class" %in% colnames(mc)) as.character(mc$class) else rep("unknown", length(te))
fam  <- if ("family" %in% colnames(mc)) as.character(mc$family) else rep("unknown", length(te))

bed7 <- tibble(
  chr = as.character(seqnames(te)),
  start = start(te) - 1L,
  end = end(te),
  strand = as.character(strand(te)),
  name = name,
  class = cls,
  family = fam
) |>
  mutate(
    chr = sub("^chr", "", chr),
    start = pmax(start, 0L),
    strand = if_else(strand %in% c("+", "-"), strand, "."),
    name = if_else(is.na(name) | name == "", paste0("TE_", row_number()), name),
    class = if_else(is.na(class) | class == "", "unknown", class),
    family = if_else(is.na(family) | family == "", "unknown", family)
  ) |>
  filter(chr %in% c(as.character(1:19), "X", "Y"), end > start) |>
  distinct() |>
  arrange(chr, start, end)

# ------------------------------------------------------------------------------
# 3. Export
# ------------------------------------------------------------------------------
out_file <- file.path(lions_ref_dir, "mm10_TE.LIONS.7col")
write_tsv(bed7, out_file, col_names = FALSE)

message("LIONS 7-col TE reference generated: ", out_file)