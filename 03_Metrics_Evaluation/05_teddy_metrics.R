# ==============================================================================
# Script: 05_teddy_metrics.R
# Purpose: Load TEDDY chimeric GTF results and compute Gene-level TP/FP/FN metrics.
# ==============================================================================


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(GenomicRanges)
})

truth_dir <- "./results"
work_dir <- "./results/TEDDY_official_by_depth"
depths <- c("5x", "10x", "25x", "50x", "100x")

# 1. Load Truth Tables
tx_truth <- read.delim(file.path(truth_dir, "official_simulated_transcript_truth_status.tsv"))
gene_truth <- read.delim(file.path(truth_dir, "official_simulated_gene_truth_status.tsv"))
iso_truth <- read.delim(file.path(truth_dir, "official_simulated_1000.isoforms.results"), check.names = FALSE)

tx_truth_expr <- tx_truth %>%
  left_join(iso_truth %>% select(transcript_id, TPM, expected_count), by = "transcript_id") %>%
  mutate(TPM = coalesce(TPM, 0), expected_count = coalesce(expected_count, 0), is_expressed_truth = TPM > 1)

truth_gene_expr <- tx_truth_expr %>% filter(is_TE_chimeric_tx, is_expressed_truth) %>% distinct(gene_name)
truth_gene_raw <- gene_truth %>% filter(is_TE_chimeric_gene) %>% distinct(gene_name)

truth_gene_expr_detail <- tx_truth_expr %>%
  filter(is_TE_chimeric_tx, is_expressed_truth) %>%
  group_by(gene_name) %>%
  summarise(
    truth_TE_tx = n_distinct(transcript_id),
    max_TPM = max(TPM, na.rm = TRUE), max_expected_count = max(expected_count, na.rm = TRUE),
    median_TPM = median(TPM, na.rm = TRUE), median_expected_count = median(expected_count, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Evaluation Loop
eval_one_depth <- function(depth) {
  meta_dir <- file.path(work_dir, depth, "meta")
  eval_dir <- file.path(work_dir, depth, "evaluation")
  dir.create(eval_dir, recursive = TRUE, showWarnings = FALSE)
  
  chi_GTF_path <- file.path(meta_dir, paste0("official_chi_GTF_", depth, ".rds"))
  if(!file.exists(chi_GTF_path)) return(NULL)
  chi_GTF <- readRDS(chi_GTF_path)
  
  pred_gene <- as.data.frame(chi_GTF) %>% transmute(gene_name = as.character(gene_name)) %>% filter(!is.na(gene_name), gene_name != "") %>% distinct()
  
  tp_gene <- inner_join(pred_gene, truth_gene_expr, by = "gene_name")
  fp_gene <- anti_join(pred_gene, truth_gene_expr, by = "gene_name")
  fn_gene <- anti_join(truth_gene_expr, pred_gene, by = "gene_name")
  
  fp_gene_detail <- fp_gene %>% mutate(predicted_by_TEDDY = TRUE) %>% left_join(truth_gene_raw %>% mutate(in_raw_TE_truth = TRUE), by = "gene_name") %>% mutate(in_raw_TE_truth = coalesce(in_raw_TE_truth, FALSE))
  fn_gene_detail <- fn_gene %>% left_join(truth_gene_expr_detail, by = "gene_name") %>% arrange(desc(max_TPM))
  tp_gene_detail <- tp_gene %>% left_join(truth_gene_expr_detail, by = "gene_name") %>% arrange(desc(max_TPM))
  
  write.table(tp_gene_detail, file.path(eval_dir, paste0("TP_gene_detail_", depth, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(fp_gene_detail, file.path(eval_dir, paste0("FP_gene_detail_", depth, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(fn_gene_detail, file.path(eval_dir, paste0("FN_gene_detail_", depth, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  
  precision <- nrow(tp_gene) / (nrow(tp_gene) + nrow(fp_gene) + 1e-9)
  recall <- nrow(tp_gene) / (nrow(tp_gene) + nrow(fn_gene) + 1e-9)
  F1 <- 2 * nrow(tp_gene) / (2 * nrow(tp_gene) + nrow(fp_gene) + nrow(fn_gene) + 1e-9)
  
  data.frame(
    depth = depth,
    chi_GTF_gene = length(unique(as.character(chi_GTF$gene_name))),
    expressed_truth_gene = nrow(truth_gene_expr), raw_truth_TE_gene = nrow(truth_gene_raw),
    TP = nrow(tp_gene), FP = nrow(fp_gene), FN = nrow(fn_gene), precision = precision, recall = recall, F1 = F1
  )
}

run_summary <- do.call(rbind, lapply(depths, eval_one_depth)) %>%
  mutate(depth_factor = factor(depth, levels = depths)) %>% arrange(depth_factor) %>% select(-depth_factor)

write.table(run_summary, file = file.path(work_dir, "TEDDY_official_by_depth_gene_level_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(run_summary, file = file.path(work_dir, "TEDDY_official_by_depth_gene_level_metrics.rds"))

message("--- TEDDY Metrics Evaluation Completed ---")
print(run_summary)