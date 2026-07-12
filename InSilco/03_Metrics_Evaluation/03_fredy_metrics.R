# ==============================================================================
# Script: 03_fredy_metrics.R
# Purpose: 
#   1. Map FREDY chimeric transcript predictions to official reference genes.
#   2. Generate FREDY-specific Truth (>50% TE overlap).
#   3. Evaluate FREDY predictions at the gene level against the FREDY50 truth.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Paths and constants (Desensitized)
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
})

data_dir <- "./data"
truth_dir <- "./results"
outbase <- "./results/FREDY_official_by_depth"

dir.create(outbase, recursive = TRUE, showWarnings = FALSE)

depths <- c("5x", "10x", "25x", "50x", "100x")
depth_dir_tbl <- tibble(
  depth = depths,
  fredy_dir = c("5x_test", "10x", "25x", "50x", "100x") 
)

min_te_fraction <- 0.5

# ------------------------------------------------------------------------------
# 1. Map FREDY predictions to 90% reference genes
# ------------------------------------------------------------------------------
message("--- Step 1: Mapping FREDY predictions to reference genes ---")
ref <- import(file.path(truth_dir, "official_simulated_reference_90pct.gtf"))
ref_gene <- ref[ref$type == "gene"]

map_one_depth <- function(depth, fredy_dir) {
  message("  -> Mapping depth: ", depth)
  chimeric_gtf <- file.path(outbase, fredy_dir, "chimeric", "chimeric.gtf")
  
  empty_res <- tibble(
    depth = depth, fredy_dir = fredy_dir, fredy_transcript_id = character(),
    fredy_gene_id = character(), chimeric_event = character(), chimeric_exon_number = character(),
    seqnames = character(), start = integer(), end = integer(), strand = character(),
    matched_gene_id = character(), matched_gene_name = character(), matched_gene_type = character(),
    overlap_width = integer()
  )
  
  if (!file.exists(chimeric_gtf) || file.info(chimeric_gtf)$size == 0) return(empty_res)
  
  fredy <- import(chimeric_gtf)
  fredy_tx <- fredy[fredy$type == "transcript"]
  if (length(fredy_tx) == 0) return(empty_res)
  
  hits <- findOverlaps(fredy_tx, ref_gene, ignore.strand = FALSE)
  hit_tbl <- tibble(
    q = queryHits(hits),
    s = subjectHits(hits),
    overlap_width = width(pintersect(fredy_tx[queryHits(hits)], ref_gene[subjectHits(hits)], ignore.strand = FALSE))
  )
  
  tx_tbl <- tibble(
    q = seq_along(fredy_tx), depth = depth, fredy_dir = fredy_dir,
    fredy_transcript_id = as.character(fredy_tx$transcript_id), fredy_gene_id = as.character(fredy_tx$gene_id),
    chimeric_event = as.character(fredy_tx$chimeric_event), chimeric_exon_number = as.character(fredy_tx$chimeric_exon_number),
    seqnames = as.character(seqnames(fredy_tx)), start = start(fredy_tx), end = end(fredy_tx),
    strand = as.character(strand(fredy_tx))
  )
  
  if (nrow(hit_tbl) == 0) {
    return(tx_tbl |> mutate(matched_gene_id = NA_character_, matched_gene_name = NA_character_, matched_gene_type = NA_character_, overlap_width = NA_integer_) |> select(-q))
  }
  
  best_hit <- hit_tbl |>
    mutate(
      matched_gene_id = as.character(ref_gene$gene_id[s]),
      matched_gene_name = as.character(ref_gene$gene_name[s]),
      matched_gene_type = as.character(ref_gene$gene_type[s])
    ) |>
    arrange(q, desc(overlap_width)) |>
    group_by(q) |> slice_head(n = 1) |> ungroup() |>
    select(q, matched_gene_id, matched_gene_name, matched_gene_type, overlap_width)
  
  tx_tbl |> left_join(best_hit, by = "q") |> select(-q)
}

fredy_tx_gene_map <- purrr::map2_dfr(depth_dir_tbl$depth, depth_dir_tbl$fredy_dir, map_one_depth)

fredy_tx_gene_map_final <- fredy_tx_gene_map |>
  mutate(
    mapping_status = if_else(!is.na(matched_gene_name) & matched_gene_name != "", "mapped_to_reference_gene", "unmapped"),
    pred_gene = if_else(mapping_status == "mapped_to_reference_gene", matched_gene_name, NA_character_)
  ) |>
  arrange(depth, mapping_status, fredy_transcript_id)

write_tsv(fredy_tx_gene_map_final, file.path(outbase, "FREDY_chimeric_transcript_gene_mapping_by_depth.tsv"))

pred_gene_by_depth <- fredy_tx_gene_map_final |>
  filter(mapping_status == "mapped_to_reference_gene") |>
  transmute(depth, gene_name = pred_gene) |>
  filter(!is.na(gene_name), gene_name != "") |>
  distinct()

# ------------------------------------------------------------------------------
# 2. Build FREDY-compatible 50% TE-overlap truth (FREDY50)
# ------------------------------------------------------------------------------
message("--- Step 2: Building FREDY50 Truth Table (>50% overlap) ---")

truth_gtf <- import(file.path(truth_dir, "official_simulated_truth_1000_transcripts.gtf"))
sim_exon <- truth_gtf[truth_gtf$type == "exon"]

te <- readRDS(file.path(data_dir, "mm10_TE.rds"))
# 注意：如果 NCBI_check 是你的自定义函数，请确保环境里已加载该函数
if(exists("NCBI_check")) te <- NCBI_check(te, ncbi_style = FALSE)

common_seq <- intersect(seqlevels(sim_exon), seqlevels(te))
sim_exon <- keepSeqlevels(sim_exon, common_seq, pruning.mode = "coarse")
te <- keepSeqlevels(te, common_seq, pruning.mode = "coarse")

sim_exon$is_TE_overlap_FREDY50 <- FALSE
sim_exon$FREDY50_TE_name <- "none"
sim_exon$FREDY50_TE_overlap_width <- 0L
sim_exon$FREDY50_TE_fraction <- 0

hits <- findOverlaps(sim_exon, te, ignore.strand = TRUE)

if (length(hits) > 0) {
  ov_width <- width(pintersect(sim_exon[queryHits(hits)], te[subjectHits(hits)], ignore.strand = TRUE))
  te_width <- width(te[subjectHits(hits)])
  te_fraction <- ov_width / te_width
  
  hit_tbl <- tibble(
    exon_idx = queryHits(hits),
    TE_name = as.character(te$name[subjectHits(hits)]),
    overlap_width = as.integer(ov_width),
    TE_fraction = as.numeric(te_fraction)
  ) |>
    filter(TE_fraction >= min_te_fraction) |>
    group_by(exon_idx) |>
    summarise(
      FREDY50_TE_name = paste(unique(TE_name), collapse = ","),
      FREDY50_TE_overlap_width = max(overlap_width, na.rm = TRUE),
      FREDY50_TE_fraction = max(TE_fraction, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (nrow(hit_tbl) > 0) {
    sim_exon$is_TE_overlap_FREDY50[hit_tbl$exon_idx] <- TRUE
    sim_exon$FREDY50_TE_name[hit_tbl$exon_idx] <- hit_tbl$FREDY50_TE_name
    sim_exon$FREDY50_TE_overlap_width[hit_tbl$exon_idx] <- hit_tbl$FREDY50_TE_overlap_width
    sim_exon$FREDY50_TE_fraction[hit_tbl$exon_idx] <- hit_tbl$FREDY50_TE_fraction
  }
}

exon_truth_FREDY50 <- as.data.frame(sim_exon) |>
  transmute(
    transcript_id = as.character(transcript_id), gene_id = as.character(gene_id),
    gene_name = as.character(gene_name), gene_type = as.character(gene_type),
    seqnames = as.character(seqnames), start, end, strand = as.character(strand),
    is_TE_overlap_FREDY50 = as.logical(is_TE_overlap_FREDY50),
    FREDY50_TE_name = as.character(FREDY50_TE_name),
    FREDY50_TE_overlap_width = as.integer(FREDY50_TE_overlap_width),
    FREDY50_TE_fraction = as.numeric(FREDY50_TE_fraction)
  )

tx_truth_FREDY50 <- exon_truth_FREDY50 |>
  group_by(transcript_id, gene_id, gene_name, gene_type) |>
  summarise(
    n_exons = n(),
    n_TE_overlap_exons_FREDY50 = sum(is_TE_overlap_FREDY50),
    is_TE_chimeric_tx_FREDY50 = any(is_TE_overlap_FREDY50),
    max_TE_overlap_width_FREDY50 = max(FREDY50_TE_overlap_width, na.rm = TRUE),
    max_TE_fraction_FREDY50 = max(FREDY50_TE_fraction, na.rm = TRUE),
    TE_names_FREDY50 = paste(unique(FREDY50_TE_name[FREDY50_TE_name != "none"]), collapse = ","),
    .groups = "drop"
  ) |> mutate(TE_names_FREDY50 = if_else(TE_names_FREDY50 == "", "none", TE_names_FREDY50))

gene_truth_FREDY50 <- tx_truth_FREDY50 |>
  group_by(gene_id, gene_name, gene_type) |>
  summarise(
    n_simulated_tx = n_distinct(transcript_id),
    n_TE_chimeric_tx_FREDY50 = sum(is_TE_chimeric_tx_FREDY50),
    is_TE_chimeric_gene_FREDY50 = any(is_TE_chimeric_tx_FREDY50),
    max_TE_fraction_FREDY50 = max(max_TE_fraction_FREDY50, na.rm = TRUE),
    .groups = "drop"
  )

iso_truth <- read.delim(file.path(truth_dir, "official_simulated_1000.isoforms.results"), check.names = FALSE)

tx_truth_expr_FREDY50 <- tx_truth_FREDY50 |>
  left_join(iso_truth |> select(transcript_id, TPM, expected_count), by = "transcript_id") |>
  mutate(
    TPM = coalesce(TPM, 0), expected_count = coalesce(expected_count, 0), is_expressed_truth = TPM > 1
  )

truth_gene_expr_FREDY50 <- tx_truth_expr_FREDY50 |> filter(is_TE_chimeric_tx_FREDY50, is_expressed_truth) |> distinct(gene_name)
truth_gene_raw_FREDY50 <- gene_truth_FREDY50 |> filter(is_TE_chimeric_gene_FREDY50) |> distinct(gene_name)

truth_gene_expr_detail_FREDY50 <- tx_truth_expr_FREDY50 |>
  filter(is_TE_chimeric_tx_FREDY50, is_expressed_truth) |>
  group_by(gene_name) |>
  summarise(
    truth_TE_tx_FREDY50 = n_distinct(transcript_id),
    max_TPM = max(TPM, na.rm = TRUE),
    max_expected_count = max(expected_count, na.rm = TRUE),
    median_TPM = median(TPM, na.rm = TRUE),
    median_expected_count = median(expected_count, na.rm = TRUE),
    max_TE_fraction_FREDY50 = max(max_TE_fraction_FREDY50, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(tx_truth_FREDY50, file.path(outbase, "FREDY50_truth_transcript_status.tsv"))
write_tsv(gene_truth_FREDY50, file.path(outbase, "FREDY50_truth_gene_status.tsv"))

# ------------------------------------------------------------------------------
# 3. Evaluate against FREDY50 truth
# ------------------------------------------------------------------------------
message("--- Step 3: Evaluating FREDY metrics ---")

eval_one_depth <- function(depth_i) {
  eval_dir <- file.path(outbase, depth_i, "evaluation_FREDY50_truth")
  dir.create(eval_dir, recursive = TRUE, showWarnings = FALSE)
  
  pred_gene <- pred_gene_by_depth |> filter(depth == depth_i) |> distinct(gene_name)
  
  tp_gene <- inner_join(pred_gene, truth_gene_expr_FREDY50, by = "gene_name")
  fp_gene <- anti_join(pred_gene, truth_gene_expr_FREDY50, by = "gene_name")
  fn_gene <- anti_join(truth_gene_expr_FREDY50, pred_gene, by = "gene_name")
  
  pred_detail <- fredy_tx_gene_map_final |>
    filter(depth == depth_i, mapping_status == "mapped_to_reference_gene") |>
    transmute(gene_name = pred_gene, fredy_transcript_id, fredy_gene_id, chimeric_event, chimeric_exon_number, overlap_width) |>
    distinct()
  
  tp_gene_detail <- tp_gene |> left_join(truth_gene_expr_detail_FREDY50, by = "gene_name") |> left_join(pred_detail, by = "gene_name") |> arrange(desc(max_TPM))
  
  fp_gene_detail <- fp_gene |>
    left_join(truth_gene_raw_FREDY50 |> mutate(in_raw_TE_truth_FREDY50 = TRUE), by = "gene_name") |>
    mutate(in_raw_TE_truth_FREDY50 = coalesce(in_raw_TE_truth_FREDY50, FALSE)) |>
    left_join(pred_detail, by = "gene_name") |> arrange(desc(in_raw_TE_truth_FREDY50), gene_name)
  
  fn_gene_detail <- fn_gene |> left_join(truth_gene_expr_detail_FREDY50, by = "gene_name") |> arrange(desc(max_TPM))
  
  write_tsv(tp_gene_detail, file.path(eval_dir, paste0("TP_gene_detail_", depth_i, ".tsv")))
  write_tsv(fp_gene_detail, file.path(eval_dir, paste0("FP_gene_detail_", depth_i, ".tsv")))
  write_tsv(fn_gene_detail, file.path(eval_dir, paste0("FN_gene_detail_", depth_i, ".tsv")))
  
  precision <- nrow(tp_gene) / (nrow(tp_gene) + nrow(fp_gene) + 1e-9)
  recall <- nrow(tp_gene) / (nrow(tp_gene) + nrow(fn_gene) + 1e-9)
  F1 <- 2 * nrow(tp_gene) / (2 * nrow(tp_gene) + nrow(fp_gene) + nrow(fn_gene) + 1e-9)
  
  data.frame(
    depth = depth_i, pred_gene = nrow(pred_gene),
    expressed_truth_gene_FREDY50 = nrow(truth_gene_expr_FREDY50), raw_truth_TE_gene_FREDY50 = nrow(truth_gene_raw_FREDY50),
    TP = nrow(tp_gene), FP = nrow(fp_gene), FN = nrow(fn_gene),
    precision = precision, recall = recall, F1 = F1,
    FP_in_raw_TE_truth_FREDY50 = sum(fp_gene_detail$in_raw_TE_truth_FREDY50),
    FP_not_in_raw_TE_truth_FREDY50 = sum(!fp_gene_detail$in_raw_TE_truth_FREDY50),
    FN_median_max_TPM = ifelse(nrow(fn_gene_detail) > 0, median(fn_gene_detail$max_TPM, na.rm = TRUE), NA_real_),
    FN_median_max_expected_count = ifelse(nrow(fn_gene_detail) > 0, median(fn_gene_detail$max_expected_count, na.rm = TRUE), NA_real_),
    FN_max_TPM = ifelse(nrow(fn_gene_detail) > 0, max(fn_gene_detail$max_TPM, na.rm = TRUE), NA_real_),
    FN_max_expected_count = ifelse(nrow(fn_gene_detail) > 0, max(fn_gene_detail$max_expected_count, na.rm = TRUE), NA_real_)
  )
}

fredy_metrics_FREDY50 <- do.call(rbind, lapply(depths, eval_one_depth)) |>
  mutate(depth_factor = factor(depth, levels = depths)) |>
  arrange(depth_factor) |> select(-depth_factor)

saveRDS(fredy_metrics_FREDY50, file.path(outbase, "FREDY_official_by_depth_gene_level_metrics_FREDY50_truth.rds"))

message("--- FREDY Evaluation Completed ---")
cat("\nFREDY50 truth summary:\n")
print(data.frame(
  raw_truth_TE_gene_FREDY50 = nrow(truth_gene_raw_FREDY50),
  expressed_truth_gene_FREDY50 = nrow(truth_gene_expr_FREDY50),
  raw_truth_TE_tx_FREDY50 = sum(tx_truth_FREDY50$is_TE_chimeric_tx_FREDY50),
  expressed_truth_TE_tx_FREDY50 = sum(tx_truth_expr_FREDY50$is_TE_chimeric_tx_FREDY50 & tx_truth_expr_FREDY50$is_expressed_truth)
))
print(fredy_metrics_FREDY50)