# ==============================================================================
# Prepare evaluable TE-host breakpoint anchors from the simulated truth GTF.
#
# Output:
#   output/simulated_truth_anchors.tsv
#   output/simulated_truth_TE_full_regions_for_TEabundance.bed
#   output/simulated_truth_TE_full_regions_for_TEabundance.tsv
# ==============================================================================


suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)

script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg)))
} else {
  normalizePath(getwd())
}

project_dir <- normalizePath(file.path(script_dir, ".."))
input_dir <- file.path(project_dir, "01_Data_Simulation", "input")
output_dir <- file.path(script_dir, "output")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gtf_file <- file.path(
  input_dir,
  "official_simulated_truth_1000_transcripts.gtf"
)

te_file <- file.path(
  input_dir,
  "mm10_TE_annotation.rds"
)

MIN_TE_OVERLAP <- 15L
SHORT_TE_THRESHOLD <- 125L
ANCHOR_WIDTH <- 50L
MIN_HOST_ANCHOR <- 15L

stopifnot(file.exists(gtf_file))
stopifnot(file.exists(te_file))

gtf <- import(gtf_file)
exons <- gtf[gtf$type == "exon"]

te <- readRDS(te_file)

if (!"name" %in% names(mcols(te)) &&
    "names" %in% names(mcols(te))) {
  te$name <- te$names
}

if (!"class" %in% names(mcols(te))) {
  te$class <- NA_character_
}

common_seq <- intersect(seqlevels(exons), seqlevels(te))

exons <- keepSeqlevels(
  exons,
  common_seq,
  pruning.mode = "coarse"
)

te <- keepSeqlevels(
  te,
  common_seq,
  pruning.mode = "coarse"
)

fmt_region <- function(chrom, start, end) {
  start <- pmax(1L, as.integer(start))
  end <- as.integer(end)
  
  ifelse(
    is.na(start) | is.na(end) | start > end,
    NA_character_,
    paste0(chrom, "|", start, "-", end)
  )
}

message("Identifying intra-exonic TE-host breakpoint anchors...")

hits <- findOverlaps(
  exons,
  te,
  minoverlap = MIN_TE_OVERLAP,
  ignore.strand = FALSE
)

inters <- pintersect(
  exons[queryHits(hits)],
  te[subjectHits(hits)],
  ignore.strand = FALSE
)

intra_raw <- tibble(
  transcript_id = as.character(
    exons[queryHits(hits)]$transcript_id
  ),
  gene_id = as.character(
    exons[queryHits(hits)]$gene_id
  ),
  gene_name = as.character(
    exons[queryHits(hits)]$gene_name
  ),
  
  TE_name = as.character(
    te[subjectHits(hits)]$name
  ),
  TE_class = as.character(
    te[subjectHits(hits)]$class
  ),
  TE_full_seqnames = as.character(
    seqnames(te[subjectHits(hits)])
  ),
  TE_full_start = start(te[subjectHits(hits)]),
  TE_full_end = end(te[subjectHits(hits)]),
  TE_full_width = width(te[subjectHits(hits)]),
  TE_full_strand = as.character(
    strand(te[subjectHits(hits)])
  ),
  
  chrom = as.character(seqnames(inters)),
  strand = as.character(
    strand(exons[queryHits(hits)])
  ),
  exon_start = start(exons[queryHits(hits)]),
  exon_end = end(exons[queryHits(hits)]),
  overlap_start = start(inters),
  overlap_end = end(inters),
  overlap_width = width(inters)
) |>
  filter(overlap_width >= MIN_TE_OVERLAP) |>
  mutate(
    left_host_width = overlap_start - exon_start,
    right_host_width = exon_end - overlap_end,
    TE_full_region = fmt_region(
      TE_full_seqnames,
      TE_full_start,
      TE_full_end
    )
  )

short_te_df <- intra_raw |>
  filter(
    overlap_width < SHORT_TE_THRESHOLD,
    left_host_width >= MIN_HOST_ANCHOR,
    right_host_width >= MIN_HOST_ANCHOR
  ) |>
  mutate(
    evidence_type = "short_full_span",
    anchor_start = overlap_start,
    anchor_end = overlap_end,
    TE_region = fmt_region(
      chrom,
      overlap_start,
      overlap_end
    ),
    Host_region = paste(
      fmt_region(
        chrom,
        pmax(exon_start, overlap_start - ANCHOR_WIDTH),
        overlap_start - 1L
      ),
      fmt_region(
        chrom,
        overlap_end + 1L,
        pmin(exon_end, overlap_end + ANCHOR_WIDTH)
      ),
      sep = ";"
    )
  )

long_te_raw <- intra_raw |>
  filter(overlap_width >= SHORT_TE_THRESHOLD)

long_left_df <- long_te_raw |>
  filter(
    overlap_start > exon_start,
    left_host_width >= MIN_HOST_ANCHOR
  ) |>
  mutate(
    evidence_type = "long_boundary_left",
    anchor_start = overlap_start,
    anchor_end = overlap_start,
    TE_region = fmt_region(
      chrom,
      overlap_start,
      pmin(
        overlap_end,
        overlap_start + ANCHOR_WIDTH - 1L
      )
    ),
    Host_region = fmt_region(
      chrom,
      pmax(
        exon_start,
        overlap_start - ANCHOR_WIDTH
      ),
      overlap_start - 1L
    )
  )

long_right_df <- long_te_raw |>
  filter(
    overlap_end < exon_end,
    right_host_width >= MIN_HOST_ANCHOR
  ) |>
  mutate(
    evidence_type = "long_boundary_right",
    anchor_start = overlap_end,
    anchor_end = overlap_end,
    TE_region = fmt_region(
      chrom,
      pmax(
        overlap_start,
        overlap_end - ANCHOR_WIDTH + 1L
      ),
      overlap_end
    ),
    Host_region = fmt_region(
      chrom,
      overlap_end + 1L,
      pmin(
        exon_end,
        overlap_end + ANCHOR_WIDTH
      )
    )
  )

all_anchor_raw <- bind_rows(
  short_te_df,
  long_left_df,
  long_right_df
) |>
  filter(
    !is.na(TE_region),
    !is.na(Host_region)
  ) |>
  arrange(
    chrom,
    anchor_start,
    anchor_end,
    transcript_id
  ) |>
  mutate(
    bp_id = paste0("BP_", row_number())
  )

all_anchors <- all_anchor_raw |>
  select(
    bp_id,
    transcript_id,
    gene_id,
    gene_name,
    TE_name,
    chrom,
    anchor_start,
    anchor_end,
    strand,
    evidence_type,
    overlap_width,
    TE_class,
    TE_region,
    Host_region
  )

write_tsv(
  all_anchors,
  file.path(
    output_dir,
    "simulated_truth_anchors.tsv"
  )
)

te_full_bed <- all_anchor_raw |>
  transmute(
    chrom = TE_full_seqnames,
    chromStart = TE_full_start - 1L,
    chromEnd = TE_full_end,
    name = paste(
      bp_id,
      transcript_id,
      gene_name,
      TE_name,
      evidence_type,
      sep = "|"
    ),
    score = 0,
    strand = TE_full_strand,
    bp_id,
    transcript_id,
    gene_id,
    gene_name,
    TE_name,
    TE_class,
    TE_full_width,
    TE_full_region,
    TE_overlap_region = TE_region,
    Host_region,
    evidence_type
  ) |>
  distinct()

write_tsv(
  te_full_bed,
  file.path(
    output_dir,
    "simulated_truth_TE_full_regions_for_TEabundance.bed"
  ),
  col_names = FALSE
)

write_tsv(
  te_full_bed,
  file.path(
    output_dir,
    "simulated_truth_TE_full_regions_for_TEabundance.tsv"
  )
)

message("Saved breakpoint-anchor files to: ", output_dir)
print(count(all_anchors, evidence_type))