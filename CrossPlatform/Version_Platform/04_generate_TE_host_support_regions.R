#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Teddy)
  library(dplyr)
  library(readr)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(rtracklayer)
})

usage <- function() {
  cat(
    paste0(
      "Usage:\n",
      "  Rscript 02_generate_TE_host_support_regions.R \\\n",
      "    <candidate_table.tsv> <combineSE.rds> <TE_annotation.rds|bed|gtf> <out_dir> \\\n",
      "    [min_TE_overlap=15] [min_host_width=15]\n\n",
      "Required candidate-table columns: transcript_id, gene_name, TE_name, TE_class\n",
      "Optional coordinate columns retained for audit: seqnames, start, end, strand\n"
    )
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) {
  usage()
  quit(status = 1L)
}

candidate_file <- args[[1]]
combine_se_file <- args[[2]]
te_file <- args[[3]]
out_dir <- args[[4]]
min_te_overlap <- if (length(args) >= 5L) as.integer(args[[5]]) else 15L
min_host_width <- if (length(args) >= 6L) as.integer(args[[6]]) else 15L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(candidate_file), file.exists(combine_se_file), file.exists(te_file))

pick_col <- function(df, candidates, required = TRUE) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0L) {
    if (required) {
      stop("Cannot find any of these columns: ", paste(candidates, collapse = ", "))
    }
    return(NA_character_)
  }
  hit[[1]]
}

read_granges <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    x <- readRDS(path)
    if (inherits(x, "GRanges")) return(x)
    if (is.data.frame(x)) return(GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = TRUE))
    stop("Unsupported RDS object: expected GRanges or data.frame")
  }
  rtracklayer::import(path)
}

candidate_raw <- readr::read_tsv(candidate_file, show_col_types = FALSE)
candidate_tx_col <- pick_col(candidate_raw, c("transcript_id", "Transcript_id", "TEDDY_vM7_tx_id"))
candidate_gene_col <- pick_col(candidate_raw, c("gene_name", "Gene_name"))
candidate_te_name_col <- pick_col(candidate_raw, c("TE_name", "TEName"))
candidate_te_class_col <- pick_col(candidate_raw, c("TE_class", "TEClass"))

candidate_tbl <- candidate_raw |>
  dplyr::mutate(
    transcript_id = as.character(.data[[candidate_tx_col]]),
    gene_name = as.character(.data[[candidate_gene_col]]),
    TE_name = as.character(.data[[candidate_te_name_col]]),
    TE_class = as.character(.data[[candidate_te_class_col]])
  )

candidate_tx <- unique(candidate_tbl$transcript_id)
combine_se <- readRDS(combine_se_file)
full_exon_gr <- Teddy::extractGTF(combineSE = combine_se, type = "exon")
te_gr <- read_granges(te_file)

full_exon_raw <- as.data.frame(full_exon_gr)
tx_col <- pick_col(full_exon_raw, c("transcript_id", "transcriptid", "tx_id"))
gene_id_col <- pick_col(full_exon_raw, c("gene_id", "geneid"))
gene_name_col <- pick_col(full_exon_raw, c("gene_name", "genename", "ref_gene_name"), required = FALSE)

full_exon_ranked <- full_exon_raw |>
  dplyr::transmute(
    seqnames = as.character(seqnames),
    start = as.integer(start),
    end = as.integer(end),
    strand = as.character(strand),
    transcript_id = as.character(.data[[tx_col]]),
    gene_id = as.character(.data[[gene_id_col]]),
    gene_name_gtf = if (!is.na(gene_name_col)) as.character(.data[[gene_name_col]]) else NA_character_
  ) |>
  dplyr::filter(
    transcript_id %in% candidate_tx,
    !is.na(seqnames), !is.na(start), !is.na(end), !is.na(transcript_id),
    start <= end
  ) |>
  dplyr::distinct(transcript_id, seqnames, start, end, strand, .keep_all = TRUE) |>
  dplyr::left_join(
    candidate_tbl |>
      dplyr::select(transcript_id, gene_name_candidate = gene_name) |>
      dplyr::distinct(),
    by = "transcript_id"
  ) |>
  dplyr::mutate(gene_name = dplyr::coalesce(gene_name_gtf, gene_name_candidate)) |>
  dplyr::group_by(transcript_id) |>
  dplyr::mutate(
    tx_strand = dplyr::case_when(
      any(strand == "-") ~ "-",
      any(strand == "+") ~ "+",
      TRUE ~ "."
    ),
    rank_key = dplyr::if_else(tx_strand == "-", -start, start),
    rank_end_key = dplyr::if_else(tx_strand == "-", -end, end)
  ) |>
  dplyr::arrange(rank_key, rank_end_key, .by_group = TRUE) |>
  dplyr::mutate(tx_exon_rank = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::select(
    transcript_id, gene_id, gene_name, seqnames, start, end, strand, tx_exon_rank
  )

exon_count_tbl <- full_exon_ranked |>
  dplyr::count(transcript_id, gene_id, gene_name, name = "n_exons_total") |>
  dplyr::mutate(eval_type = dplyr::if_else(n_exons_total >= 2L, "multi_exon", "single_exon"))

first_exon_tbl <- full_exon_ranked |>
  dplyr::filter(tx_exon_rank == 1L) |>
  dplyr::left_join(
    exon_count_tbl |>
      dplyr::select(transcript_id, n_exons_total, eval_type),
    by = "transcript_id"
  ) |>
  dplyr::left_join(
    candidate_tbl |>
      dplyr::select(
        transcript_id,
        candidate_TE_name = TE_name,
        candidate_TE_class = TE_class
      ) |>
      dplyr::distinct(),
    by = "transcript_id"
  )

second_exon_tbl <- full_exon_ranked |>
  dplyr::filter(tx_exon_rank == 2L) |>
  dplyr::transmute(
    transcript_id,
    exon2_seqnames = seqnames,
    exon2_start = start,
    exon2_end = end,
    exon2_strand = strand
  )

first_exon_gr <- GenomicRanges::GRanges(
  seqnames = first_exon_tbl$seqnames,
  ranges = IRanges::IRanges(first_exon_tbl$start, first_exon_tbl$end),
  strand = first_exon_tbl$strand
)
S4Vectors::mcols(first_exon_gr) <- S4Vectors::DataFrame(
  transcript_id = first_exon_tbl$transcript_id,
  gene_id = first_exon_tbl$gene_id,
  gene_name = first_exon_tbl$gene_name,
  n_exons_total = first_exon_tbl$n_exons_total,
  eval_type = first_exon_tbl$eval_type,
  candidate_TE_name = first_exon_tbl$candidate_TE_name,
  candidate_TE_class = first_exon_tbl$candidate_TE_class
)

te_raw <- as.data.frame(te_gr)
te_name_col <- pick_col(te_raw, c("names", "name", "repName", "TE_name"), required = FALSE)
te_class_col <- pick_col(te_raw, c("class", "repClass", "TE_class"), required = FALSE)

ov <- GenomicRanges::findOverlaps(
  first_exon_gr,
  te_gr,
  ignore.strand = TRUE,
  minoverlap = min_te_overlap
)

if (length(ov) == 0L) stop("No TE overlaps were found for candidate first exons")

first_hit <- first_exon_gr[S4Vectors::queryHits(ov)]
te_hit <- te_gr[S4Vectors::subjectHits(ov)]
te_overlap_gr <- GenomicRanges::pintersect(first_hit, te_hit, ignore.strand = TRUE)

te_pair_tbl <- tibble::tibble(
  object_index = S4Vectors::queryHits(ov),
  te_index = S4Vectors::subjectHits(ov),
  transcript_id = as.character(first_hit$transcript_id),
  gene_id = as.character(first_hit$gene_id),
  gene_name = as.character(first_hit$gene_name),
  eval_type = as.character(first_hit$eval_type),
  n_exons_total = as.integer(first_hit$n_exons_total),
  seqnames = as.character(GenomicRanges::seqnames(first_hit)),
  strand = as.character(GenomicRanges::strand(first_hit)),
  exon1_start = GenomicRanges::start(first_hit),
  exon1_end = GenomicRanges::end(first_hit),
  te_start = GenomicRanges::start(te_overlap_gr),
  te_end = GenomicRanges::end(te_overlap_gr),
  TE_name = if (!is.na(te_name_col)) as.character(S4Vectors::mcols(te_hit)[[te_name_col]]) else as.character(first_hit$candidate_TE_name),
  TE_class = if (!is.na(te_class_col)) as.character(S4Vectors::mcols(te_hit)[[te_class_col]]) else as.character(first_hit$candidate_TE_class)
) |>
  dplyr::mutate(
    te_len = te_end - te_start + 1L,
    left_host_start = exon1_start,
    left_host_end = te_start - 1L,
    left_host_len = left_host_end - left_host_start + 1L,
    right_host_start = te_end + 1L,
    right_host_end = exon1_end,
    right_host_len = right_host_end - right_host_start + 1L
  ) |>
  dplyr::left_join(second_exon_tbl, by = "transcript_id")

same_left_pair <- te_pair_tbl |>
  dplyr::filter(left_host_len >= min_host_width) |>
  dplyr::transmute(
    pair_id = paste(transcript_id, object_index, te_index, "same_exon_left", sep = "__"),
    transcript_id, gene_id, gene_name, eval_type, n_exons_total, TE_name, TE_class,
    te_seqnames = seqnames, te_start, te_end, te_strand = strand,
    host_seqnames = seqnames, host_start = left_host_start, host_end = left_host_end,
    host_strand = strand, host_anchor_type = "same_exon_left_flank"
  )

same_right_pair <- te_pair_tbl |>
  dplyr::filter(right_host_len >= min_host_width) |>
  dplyr::transmute(
    pair_id = paste(transcript_id, object_index, te_index, "same_exon_right", sep = "__"),
    transcript_id, gene_id, gene_name, eval_type, n_exons_total, TE_name, TE_class,
    te_seqnames = seqnames, te_start, te_end, te_strand = strand,
    host_seqnames = seqnames, host_start = right_host_start, host_end = right_host_end,
    host_strand = strand, host_anchor_type = "same_exon_right_flank"
  )

exon2_pair <- te_pair_tbl |>
  dplyr::filter(
    eval_type == "multi_exon",
    !is.na(exon2_seqnames), !is.na(exon2_start), !is.na(exon2_end)
  ) |>
  dplyr::transmute(
    pair_id = paste(transcript_id, object_index, te_index, "exon2", sep = "__"),
    transcript_id, gene_id, gene_name, eval_type, n_exons_total, TE_name, TE_class,
    te_seqnames = seqnames, te_start, te_end, te_strand = strand,
    host_seqnames = exon2_seqnames, host_start = exon2_start, host_end = exon2_end,
    host_strand = exon2_strand, host_anchor_type = "downstream_exon2"
  )

pair_meta <- dplyr::bind_rows(same_left_pair, same_right_pair, exon2_pair) |>
  dplyr::mutate(
    te_len = te_end - te_start + 1L,
    host_len = host_end - host_start + 1L,
    pair_len = te_len + host_len
  ) |>
  dplyr::filter(te_len >= min_te_overlap, host_len >= min_host_width) |>
  dplyr::distinct() |>
  dplyr::select(
    pair_id, transcript_id, gene_id, gene_name, eval_type, n_exons_total,
    TE_name, TE_class,
    te_seqnames, te_start, te_end, te_strand,
    host_seqnames, host_start, host_end, host_strand,
    host_anchor_type, te_len, host_len, pair_len
  )

te_anchor_bed <- pair_meta |>
  dplyr::transmute(
    chrom = te_seqnames,
    chromStart = pmax(te_start - 1L, 0L),
    chromEnd = te_end,
    name = paste(pair_id, transcript_id, gene_name, "TE_anchor", host_anchor_type, sep = "|"),
    score = 0L,
    strand = dplyr::if_else(te_strand %in% c("+", "-"), te_strand, ".")
  )

host_anchor_bed <- pair_meta |>
  dplyr::transmute(
    chrom = host_seqnames,
    chromStart = pmax(host_start - 1L, 0L),
    chromEnd = host_end,
    name = paste(pair_id, transcript_id, gene_name, "host_anchor", host_anchor_type, sep = "|"),
    score = 0L,
    strand = dplyr::if_else(host_strand %in% c("+", "-"), host_strand, ".")
  )

anchor_bed <- dplyr::bind_rows(te_anchor_bed, host_anchor_bed) |>
  dplyr::distinct()

evaluable_tx <- unique(pair_meta$transcript_id)
unevaluable_tbl <- candidate_tbl |>
  dplyr::filter(!transcript_id %in% evaluable_tx)

pair_meta_file <- file.path(out_dir, "TE_host_pair_meta.tsv")
anchor_bed_file <- file.path(out_dir, "TE_host_eval_anchors.bed")
unevaluable_file <- file.path(out_dir, "TE_host_unevaluable_candidates.tsv")
summary_file <- file.path(out_dir, "TE_host_region_generation_summary.tsv")

readr::write_tsv(pair_meta, pair_meta_file)
readr::write_tsv(anchor_bed, anchor_bed_file, col_names = FALSE)
readr::write_tsv(unevaluable_tbl, unevaluable_file)

summary_tbl <- tibble::tibble(
  metric = c(
    "n_candidates",
    "n_multi_exon_candidates",
    "n_single_exon_candidates",
    "n_evaluable_candidates",
    "n_unevaluable_candidates",
    "n_TE_host_pairs"
  ),
  value = c(
    length(candidate_tx),
    sum(exon_count_tbl$eval_type == "multi_exon"),
    sum(exon_count_tbl$eval_type == "single_exon"),
    length(evaluable_tx),
    nrow(unevaluable_tbl),
    nrow(pair_meta)
  )
)
readr::write_tsv(summary_tbl, summary_file)

cat("Generated:\n")
cat("  ", pair_meta_file, "\n", sep = "")
cat("  ", anchor_bed_file, "\n", sep = "")
cat("  ", unevaluable_file, "\n", sep = "")
cat("  ", summary_file, "\n", sep = "")
print(summary_tbl)
