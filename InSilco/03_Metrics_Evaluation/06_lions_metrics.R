# ==============================================================================
# Script: 06_lions_metrics.R
# Purpose: Extract final LIONS predictions from .lion files and evaluate accuracy against
#          a strict TE-initiated ground truth.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

truth_dir <- "./results"
outbase <- file.path(truth_dir, "LIONS_official_by_depth")
depths <- c("5x", "10x", "25x", "50x", "100x")

# ------------------------------------------------------------------------------
# 1. Truth Tables
# ------------------------------------------------------------------------------
tx_truth <- read.delim(file.path(truth_dir, "official_simulated_transcript_truth_status.tsv"))
exon_truth <- read_tsv(file.path(truth_dir, "official_simulated_exon_truth_status.tsv"), show_col_types = FALSE)
iso_truth <- read.delim(file.path(truth_dir, "official_simulated_1000.isoforms.results"), check.names = FALSE)

# ------------------------------------------------------------------------------
# 2. Build TE-initiated-compatible truth
#    Rule: 5' exon of transcript must be TE-overlap exon, and TPM > 1
# ------------------------------------------------------------------------------
TEinitiated_tx_truth <- exon_truth |>
  mutate(is_TE_overlap_exon = as.logical(is_TE_overlap_exon)) |>
  group_by(transcript_id) |>
  filter(
    (strand == "+" & start == min(start)) |
      (strand == "-" & end == max(end))
  ) |>
  summarise(gene_name = dplyr::first(gene_name), TEinitiated_compatible_tx = any(is_TE_overlap_exon, na.rm = TRUE), .groups = "drop")

TEinitiated_tx_truth_expr <- tx_truth |>
  select(transcript_id, gene_name) |>
  distinct() |>
  left_join(TEinitiated_tx_truth, by = c("transcript_id", "gene_name")) |>
  left_join(iso_truth |> select(transcript_id, TPM, expected_count), by = "transcript_id") |>
  mutate(
    TEinitiated_compatible_tx = dplyr::coalesce(TEinitiated_compatible_tx, FALSE),
    TPM = dplyr::coalesce(TPM, 0), expected_count = dplyr::coalesce(expected_count, 0),
    is_expressed_truth = TPM > 1
  )

truth_gene_expr_TEinitiated <- TEinitiated_tx_truth_expr |> filter(TEinitiated_compatible_tx, is_expressed_truth) |> distinct(gene_name)

truth_gene_expr_TEinitiated_detail <- TEinitiated_tx_truth_expr |>
  filter(TEinitiated_compatible_tx, is_expressed_truth) |>
  group_by(gene_name) |>
  summarise(
    truth_TE_tx_TEinitiated = n_distinct(transcript_id),
    max_TPM = max(TPM, na.rm = TRUE), max_expected_count = max(expected_count, na.rm = TRUE),
    median_TPM = median(TPM, na.rm = TRUE), median_expected_count = median(expected_count, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 3. Read final chimSort/simuOptimal LIONS predictions
# ------------------------------------------------------------------------------
read_lions_pred_gene <- function(depth_i) {
  f <- file.path(outbase, depth_i, paste0("official_simulated_", depth_i, "_noise0.1.lion"))
  if (!file.exists(f) || file.info(f)$size == 0) {
    warning("Missing or empty LIONS output: ", f)
    return(tibble(depth = depth_i, gene_name = character()))
  }
  
  x <- read_tsv(f, show_col_types = FALSE, col_types = cols(.default = col_character()))
  stopifnot("transcriptID" %in% names(x))

  tibble(depth = depth_i, transcriptID = x$transcriptID) |>
    mutate(gene_name = sub(":.*$", "", transcriptID), gene_name = trimws(gene_name)) |>
    filter(!is.na(gene_name), gene_name != "", gene_name != "NA", gene_name != "None") |>
    distinct(depth, gene_name)
}

pred_gene_by_depth <- bind_rows(lapply(depths, read_lions_pred_gene))
write_tsv(pred_gene_by_depth, file.path(outbase, "LIONS_pred_gene_by_depth.tsv"))

# ------------------------------------------------------------------------------
# 4. Evaluate all depths
# ------------------------------------------------------------------------------
eval_one_depth <- function(depth_i) {
  eval_dir <- file.path(outbase, depth_i, "evaluation_gene_level_TEinitiated_truth")
  dir.create(eval_dir, recursive = TRUE, showWarnings = FALSE)
  
  pred_gene <- pred_gene_by_depth |> filter(depth == depth_i) |> distinct(gene_name)
  
  tp_gene <- inner_join(pred_gene, truth_gene_expr_TEinitiated, by = "gene_name")
  fp_gene <- anti_join(pred_gene, truth_gene_expr_TEinitiated, by = "gene_name")
  fn_gene <- anti_join(truth_gene_expr_TEinitiated, pred_gene, by = "gene_name")
  
  tp_gene_detail <- tp_gene |> left_join(truth_gene_expr_TEinitiated_detail, by = "gene_name") |> arrange(desc(max_TPM))
  fp_gene_detail <- fp_gene |> arrange(gene_name)
  fn_gene_detail <- fn_gene |> left_join(truth_gene_expr_TEinitiated_detail, by = "gene_name") |> arrange(desc(max_TPM))
  
  write_tsv(tp_gene_detail, file.path(eval_dir, paste0("TP_gene_detail_", depth_i, ".tsv")))
  write_tsv(fp_gene_detail, file.path(eval_dir, paste0("FP_gene_detail_", depth_i, ".tsv")))
  write_tsv(fn_gene_detail, file.path(eval_dir, paste0("FN_gene_detail_", depth_i, ".tsv")))
  
  tp <- nrow(tp_gene); fp <- nrow(fp_gene); fn <- nrow(fn_gene)
  precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  F1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
  
  tibble(
    depth = depth_i, pred_gene = nrow(pred_gene), expressed_truth_gene_TEinitiated = nrow(truth_gene_expr_TEinitiated),
    TP = tp, FP = fp, FN = fn, precision = precision, recall = recall, F1 = F1,
    FN_median_max_TPM = ifelse(nrow(fn_gene_detail) > 0, median(fn_gene_detail$max_TPM, na.rm = TRUE), NA_real_),
    FN_median_max_expected_count = ifelse(nrow(fn_gene_detail) > 0, median(fn_gene_detail$max_expected_count, na.rm = TRUE), NA_real_),
    FN_max_TPM = ifelse(nrow(fn_gene_detail) > 0, max(fn_gene_detail$max_TPM, na.rm = TRUE), NA_real_),
    FN_max_expected_count = ifelse(nrow(fn_gene_detail) > 0, max(fn_gene_detail$max_expected_count, na.rm = TRUE), NA_real_)
  )
}

lions_metrics_TEinitiated <- bind_rows(lapply(depths, eval_one_depth)) |>
  mutate(depth = factor(depth, levels = depths)) |> arrange(depth)

write_tsv(lions_metrics_TEinitiated, file.path(outbase, "LIONS_official_by_depth_gene_level_metrics_TEinitiated_truth.tsv"))
saveRDS(lions_metrics_TEinitiated, file.path(outbase, "LIONS_official_by_depth_gene_level_metrics_TEinitiated_truth.rds"))

message("--- LIONS Metrics Evaluation Completed ---")
cat("Truth genes, TE-initiated compatible: ", nrow(truth_gene_expr_TEinitiated), "\n")
print(lions_metrics_TEinitiated)
