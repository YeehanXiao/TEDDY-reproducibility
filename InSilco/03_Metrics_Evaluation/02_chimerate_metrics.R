#!/usr/bin/env Rscript
# ==============================================================================
# Script: 02_chimerate_metrics.R
# Purpose: Evaluate ChimeraTE predictions against simulated truth at multiple depths.
#          Calculates Precision, Recall, and F1 at the gene level.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Paths and constants (Desensitized)
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

truth_dir <- "./results"                                      
chim_dir <- "./results/ChimeraTE_official_by_depth"          
eval_dir <- file.path(chim_dir, "evaluation")                 
dir.create(eval_dir, recursive = TRUE, showWarnings = FALSE)

depths <- c("5x", "10x", "25x", "50x", "100x")
strands <- c("fwd-stranded", "rf-stranded")
chimera_types <- c("TE-exonized", "TE-initiated", "TE-terminated")

# ------------------------------------------------------------------------------
# 1. Load Truth Tables
# ------------------------------------------------------------------------------
tx_truth <- read.delim(file.path(truth_dir, "official_simulated_transcript_truth_status.tsv"))
gene_truth <- read.delim(file.path(truth_dir, "official_simulated_gene_truth_status.tsv"))
iso_truth <- read.delim(file.path(truth_dir, "official_simulated_1000.isoforms.results"), check.names = FALSE)

tx_truth_expr <- tx_truth %>%
  left_join(
    iso_truth %>% select(transcript_id, TPM, expected_count),
    by = "transcript_id"
  ) %>%
  mutate(
    TPM = coalesce(TPM, 0),
    expected_count = coalesce(expected_count, 0),
    is_expressed_truth = TPM > 1
  )

truth_gene_expr <- tx_truth_expr %>%
  filter(is_TE_chimeric_tx, is_expressed_truth) %>%
  distinct(gene_id, gene_name)

truth_gene_raw <- gene_truth %>%
  filter(is_TE_chimeric_gene) %>%
  distinct(gene_id, gene_name)

truth_gene_expr_detail <- tx_truth_expr %>%
  filter(is_TE_chimeric_tx, is_expressed_truth) %>%
  group_by(gene_id, gene_name) %>%
  summarise(
    truth_TE_tx = n_distinct(transcript_id),
    max_TPM = max(TPM, na.rm = TRUE),
    max_expected_count = max(expected_count, na.rm = TRUE),
    median_TPM = median(TPM, na.rm = TRUE),
    median_expected_count = median(expected_count, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 2. Read ChimeraTE Predictions
# ------------------------------------------------------------------------------
read_chimerate_one <- function(depth, strand, chimera_type) {
  f <- file.path(
    chim_dir,
    paste0(depth, "_", strand),
    paste0("rep1_", depth),
    paste0(chimera_type, "-rep1_", depth, ".tsv")
  )
  
  if (!file.exists(f) || file.info(f)$size == 0) {
    return(tibble())
  }
  
  x <- suppressWarnings(read_tsv(
    f,
    col_names = FALSE,
    show_col_types = FALSE,
    progress = FALSE
  ))
  
  if (nrow(x) == 0 || ncol(x) < 1) {
    return(tibble())
  }
  
  colnames(x)[seq_len(min(ncol(x), 8))] <- c(
    "gene_id", "gene_strand", "gene_coord",
    "TE_name", "TE_strand", "TE_coord",
    "support_or_overlap", "chimera_subtype"
  )[seq_len(min(ncol(x), 8))]
  
  x %>%
    mutate(
      depth = depth,
      strand_mode = strand,
      chimera_type = chimera_type,
      source_file = f,
      gene_id = as.character(gene_id)
    ) %>%
    filter(!is.na(gene_id), gene_id != "")
}

read_chimerate_depth <- function(depth) {
  bind_rows(lapply(strands, function(strand) {
    bind_rows(lapply(chimera_types, function(chimera_type) {
      read_chimerate_one(depth, strand, chimera_type)
    }))
  }))
}

# ------------------------------------------------------------------------------
# 3. Generic Evaluator
# ------------------------------------------------------------------------------
evaluate_depth <- function(depth) {
  pred_raw <- read_chimerate_depth(depth)
  
  pred_gene <- pred_raw %>%
    distinct(gene_id) %>%
    left_join(
      gene_truth %>% distinct(gene_id, gene_name),
      by = "gene_id"
    )
  
  tp_gene <- inner_join(pred_gene, truth_gene_expr, by = c("gene_id", "gene_name"))
  fp_gene <- anti_join(pred_gene, truth_gene_expr, by = c("gene_id", "gene_name"))
  fn_gene <- anti_join(truth_gene_expr, pred_gene, by = c("gene_id", "gene_name"))
  
  fp_gene_detail <- fp_gene %>%
    mutate(predicted_by_ChimeraTE = TRUE) %>%
    left_join(
      truth_gene_raw %>% mutate(in_raw_TE_truth = TRUE),
      by = c("gene_id", "gene_name")
    ) %>%
    mutate(in_raw_TE_truth = coalesce(in_raw_TE_truth, FALSE))
  
  tp_gene_detail <- tp_gene %>%
    left_join(truth_gene_expr_detail, by = c("gene_id", "gene_name")) %>%
    arrange(desc(max_TPM))
  
  fn_gene_detail <- fn_gene %>%
    left_join(truth_gene_expr_detail, by = c("gene_id", "gene_name")) %>%
    arrange(desc(max_TPM))
  
  pred_raw_detail <- pred_raw %>%
    left_join(
      gene_truth %>% distinct(gene_id, gene_name),
      by = "gene_id"
    ) %>%
    arrange(chimera_type, gene_id, TE_name)
  
  # Export detail tables
  write.table(pred_raw_detail, file.path(eval_dir, paste0("ChimeraTE_pred_raw_detail_", depth, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(tp_gene_detail, file.path(eval_dir, paste0("ChimeraTE_TP_gene_detail_", depth, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(fp_gene_detail, file.path(eval_dir, paste0("ChimeraTE_FP_gene_detail_", depth, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(fn_gene_detail, file.path(eval_dir, paste0("ChimeraTE_FN_gene_detail_", depth, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  
  TP <- nrow(tp_gene)
  FP <- nrow(fp_gene)
  FN <- nrow(fn_gene)
  
  data.frame(
    depth = depth,
    pred_event_raw = nrow(pred_raw),
    pred_gene = nrow(pred_gene),
    expressed_truth_gene = nrow(truth_gene_expr),
    raw_truth_TE_gene = nrow(truth_gene_raw),
    TP = TP,
    FP = FP,
    FN = FN,
    precision = TP / (TP + FP + 1e-9),
    recall = TP / (TP + FN + 1e-9),
    F1 = 2 * TP / (2 * TP + FP + FN + 1e-9),
    FP_in_raw_TE_truth = sum(fp_gene_detail$in_raw_TE_truth),
    FP_not_in_raw_TE_truth = sum(!fp_gene_detail$in_raw_TE_truth),
    pred_TE_exonized_events = sum(pred_raw$chimera_type == "TE-exonized"),
    pred_TE_initiated_events = sum(pred_raw$chimera_type == "TE-initiated"),
    pred_TE_terminated_events = sum(pred_raw$chimera_type == "TE-terminated")
  )
}

# ------------------------------------------------------------------------------
# 4. Evaluate and Summarize Metrics
# ------------------------------------------------------------------------------
chimerate_metrics_gene <- bind_rows(lapply(depths, evaluate_depth))

chimerate_metrics_gene <- chimerate_metrics_gene %>%
  mutate(depth = factor(depth, levels = depths)) %>%
  arrange(depth)

write.table(
  chimerate_metrics_gene,
  file.path(eval_dir, "ChimeraTE_official_by_depth_gene_level_metrics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  chimerate_metrics_gene,
  file.path(eval_dir, "ChimeraTE_official_by_depth_gene_level_metrics.rds")
)

message("--- ChimeraTE Evaluation Completed ---")
message("Metrics saved to: ", eval_dir)
