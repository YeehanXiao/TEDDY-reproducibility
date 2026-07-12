# ==============================================================================
# Script: 01_teprof2_metrics.R
# Purpose: Evaluate TEProf2 predictions against simulated truth at multiple depths.
#          Calculates Precision, Recall, and F1 at the gene level.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Paths and constants (Desensitized)
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

truth_dir <- "./results"                     
benchmark_dir <- "./results/TEProf2_benchmark" 

depths <- c("5x", "10x", "25x", "50x", "100x")

dir.create(benchmark_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Truth tables
# ------------------------------------------------------------------------------
tx_truth <- read.delim(file.path(truth_dir, "official_simulated_transcript_truth_status.tsv"))
gene_truth <- read.delim(file.path(truth_dir, "official_simulated_gene_truth_status.tsv"))
exon_truth <- read_tsv(file.path(truth_dir, "official_simulated_exon_truth_status.tsv"), show_col_types = FALSE)
iso_truth <- read.delim(file.path(truth_dir, "official_simulated_1000.isoforms.results"), check.names = FALSE)

tx_truth_expr <- tx_truth |>
  left_join(
    iso_truth |> select(transcript_id, TPM, expected_count),
    by = "transcript_id"
  ) |>
  mutate(
    TPM = dplyr::coalesce(TPM, 0),
    expected_count = dplyr::coalesce(expected_count, 0),
    is_expressed_truth = TPM > 1
  )

# main truth: any TE-chimeric transcript + TPM > 1
truth_gene_expr_main <- tx_truth_expr |>
  filter(is_TE_chimeric_tx, is_expressed_truth) |>
  distinct(gene_name)

truth_gene_raw_main <- gene_truth |>
  filter(is_TE_chimeric_gene) |>
  distinct(gene_name)

truth_gene_expr_main_detail <- tx_truth_expr |>
  filter(is_TE_chimeric_tx, is_expressed_truth) |>
  group_by(gene_name) |>
  summarise(
    truth_TE_tx = n_distinct(transcript_id),
    max_TPM = max(TPM, na.rm = TRUE),
    max_expected_count = max(expected_count, na.rm = TRUE),
    median_TPM = median(TPM, na.rm = TRUE),
    median_expected_count = median(expected_count, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 2. TEProf2-compatible truth
#    5' exon of transcript must be TE-overlap exon, then TPM > 1
# ------------------------------------------------------------------------------
teprof2_tx_truth <- exon_truth |>
  mutate(
    tx_exon_rank = as.integer(tx_exon_rank),
    is_TE_overlap_exon = as.logical(is_TE_overlap_exon)
  ) |>
  group_by(transcript_id) |>
  filter(
    (strand == "+" & tx_exon_rank == min(tx_exon_rank, na.rm = TRUE)) |
      (strand == "-" & tx_exon_rank == max(tx_exon_rank, na.rm = TRUE))
  ) |>
  summarise(
    gene_name = dplyr::first(gene_name),
    teprof2_compatible_tx = any(is_TE_overlap_exon, na.rm = TRUE),
    .groups = "drop"
  )

teprof2_tx_truth_expr <- tx_truth |>
  select(transcript_id, gene_name) |>
  distinct() |>
  left_join(teprof2_tx_truth, by = c("transcript_id", "gene_name")) |>
  left_join(
    iso_truth |> select(transcript_id, TPM, expected_count),
    by = "transcript_id"
  ) |>
  mutate(
    teprof2_compatible_tx = dplyr::coalesce(teprof2_compatible_tx, FALSE),
    TPM = dplyr::coalesce(TPM, 0),
    expected_count = dplyr::coalesce(expected_count, 0),
    is_expressed_truth = TPM > 1
  )

truth_gene_expr_teprof2 <- teprof2_tx_truth_expr |>
  filter(teprof2_compatible_tx, is_expressed_truth) |>
  distinct(gene_name)

truth_gene_expr_teprof2_detail <- teprof2_tx_truth_expr |>
  filter(teprof2_compatible_tx, is_expressed_truth) |>
  group_by(gene_name) |>
  summarise(
    truth_TE_tx_teprof2 = n_distinct(transcript_id),
    max_TPM = max(TPM, na.rm = TRUE),
    max_expected_count = max(expected_count, na.rm = TRUE),
    median_TPM = median(TPM, na.rm = TRUE),
    median_expected_count = median(expected_count, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 3. Read TEProf2 predictions
#    Format confirmed: no header; prefer X16, fallback X3
# ------------------------------------------------------------------------------
read_teprof2_pred_gene <- function(depth_i) {
  f <- file.path(
    benchmark_dir, paste0(depth_i, "_test"), "assembly",
    paste0("official_simulated_", depth_i, "_noise0.1.stringtie.gtf_annotated_filtered_test_all")
  )
  
  if (!file.exists(f) || file.info(f)$size == 0) {
    return(tibble(depth = depth_i, gene_name = character()))
  }
  
  x <- read_tsv(
    f,
    col_names = FALSE,
    show_col_types = FALSE
  )
  
  if (ncol(x) < 16) {
    stop("Unexpected TEProf2 column number in: ", f)
  }
  
  x |>
    transmute(
      depth = depth_i,
      gene_name = dplyr::coalesce(X16, X3)
    ) |>
    mutate(
      gene_name = trimws(gene_name)
    ) |>
    filter(!is.na(gene_name), gene_name != "", gene_name != "None") |>
    distinct()
}

pred_gene_by_depth <- bind_rows(lapply(depths, read_teprof2_pred_gene))

write_tsv(
  pred_gene_by_depth,
  file.path(benchmark_dir, "TEProf2_pred_gene_by_depth.tsv")
)

# ------------------------------------------------------------------------------
# 4. Generic evaluator
# ------------------------------------------------------------------------------
eval_one_depth <- function(
    depth_i,
    pred_gene_by_depth,
    truth_gene_expr,
    truth_gene_expr_detail,
    truth_gene_raw = NULL,
    eval_subdir = "evaluation_gene_level",
    truth_label = "main"
) {
  eval_dir <- file.path(benchmark_dir, paste0(depth_i, "_test"), eval_subdir)
  dir.create(eval_dir, recursive = TRUE, showWarnings = FALSE)
  
  pred_gene <- pred_gene_by_depth |>
    filter(depth == depth_i) |>
    distinct(gene_name)
  
  tp_gene <- inner_join(pred_gene, truth_gene_expr, by = "gene_name")
  fp_gene <- anti_join(pred_gene, truth_gene_expr, by = "gene_name")
  fn_gene <- anti_join(truth_gene_expr, pred_gene, by = "gene_name")
  
  tp_gene_detail <- tp_gene |>
    left_join(truth_gene_expr_detail, by = "gene_name") |>
    arrange(desc(max_TPM))
  
  fp_gene_detail <- fp_gene
  
  if (!is.null(truth_gene_raw)) {
    fp_gene_detail <- fp_gene_detail |>
      left_join(
        truth_gene_raw |> mutate(in_raw_TE_truth = TRUE),
        by = "gene_name"
      ) |>
      mutate(in_raw_TE_truth = dplyr::coalesce(in_raw_TE_truth, FALSE)) |>
      arrange(desc(in_raw_TE_truth), gene_name)
  } else {
    fp_gene_detail <- fp_gene_detail |>
      arrange(gene_name)
  }
  
  fn_gene_detail <- fn_gene |>
    left_join(truth_gene_expr_detail, by = "gene_name") |>
    arrange(desc(max_TPM))
  
  write_tsv(tp_gene_detail, file.path(eval_dir, paste0("TP_gene_detail_", depth_i, ".tsv")))
  write_tsv(fp_gene_detail, file.path(eval_dir, paste0("FP_gene_detail_", depth_i, ".tsv")))
  write_tsv(fn_gene_detail, file.path(eval_dir, paste0("FN_gene_detail_", depth_i, ".tsv")))
  
  tp <- nrow(tp_gene)
  fp <- nrow(fp_gene)
  fn <- nrow(fn_gene)
  
  precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  F1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) {
    0
  } else {
    2 * precision * recall / (precision + recall)
  }
  
  out <- tibble(
    depth = depth_i,
    pred_gene = nrow(pred_gene),
    TP = tp,
    FP = fp,
    FN = fn,
    precision = precision,
    recall = recall,
    F1 = F1
  )
  
  if (truth_label == "main") {
    out <- out |>
      mutate(
        expressed_truth_gene = nrow(truth_gene_expr),
        raw_truth_TE_gene = nrow(truth_gene_raw),
        FP_in_raw_TE_truth = sum(fp_gene_detail$in_raw_TE_truth),
        FP_not_in_raw_TE_truth = sum(!fp_gene_detail$in_raw_TE_truth),
        FN_median_max_TPM = ifelse(nrow(fn_gene_detail) > 0, median(fn_gene_detail$max_TPM, na.rm = TRUE), NA_real_),
        FN_median_max_expected_count = ifelse(nrow(fn_gene_detail) > 0, median(fn_gene_detail$max_expected_count, na.rm = TRUE), NA_real_),
        FN_max_TPM = ifelse(nrow(fn_gene_detail) > 0, max(fn_gene_detail$max_TPM, na.rm = TRUE), NA_real_),
        FN_max_expected_count = ifelse(nrow(fn_gene_detail) > 0, max(fn_gene_detail$max_expected_count, na.rm = TRUE), NA_real_)
      ) |>
      select(
        depth, pred_gene, expressed_truth_gene, raw_truth_TE_gene,
        TP, FP, FN, precision, recall, F1,
        FP_in_raw_TE_truth, FP_not_in_raw_TE_truth,
        FN_median_max_TPM, FN_median_max_expected_count,
        FN_max_TPM, FN_max_expected_count
      )
  } else if (truth_label == "teprof2") {
    out <- out |>
      mutate(
        expressed_truth_gene_TEProf2 = nrow(truth_gene_expr)
      ) |>
      select(
        depth, pred_gene, expressed_truth_gene_TEProf2,
        TP, FP, FN, precision, recall, F1
      )
  }
  
  out
}

# ------------------------------------------------------------------------------
# 5. Evaluate under main truth
# ------------------------------------------------------------------------------
teprof2_metrics_main <- bind_rows(lapply(depths, function(depth_i) {
  eval_one_depth(
    depth_i = depth_i,
    pred_gene_by_depth = pred_gene_by_depth,
    truth_gene_expr = truth_gene_expr_main,
    truth_gene_expr_detail = truth_gene_expr_main_detail,
    truth_gene_raw = truth_gene_raw_main,
    eval_subdir = "evaluation_gene_level",
    truth_label = "main"
  )
})) |>
  mutate(depth = factor(depth, levels = depths)) |>
  arrange(depth)

write_tsv(teprof2_metrics_main, file.path(benchmark_dir, "TEProf2_official_by_depth_gene_level_metrics.tsv"))
saveRDS(teprof2_metrics_main, file.path(benchmark_dir, "TEProf2_official_by_depth_gene_level_metrics.rds"))

# ------------------------------------------------------------------------------
# 6. Evaluate under TEProf2-compatible truth
# ------------------------------------------------------------------------------
teprof2_metrics_teinit <- bind_rows(lapply(depths, function(depth_i) {
  eval_one_depth(
    depth_i = depth_i,
    pred_gene_by_depth = pred_gene_by_depth,
    truth_gene_expr = truth_gene_expr_teprof2,
    truth_gene_expr_detail = truth_gene_expr_teprof2_detail,
    truth_gene_raw = NULL,
    eval_subdir = "evaluation_gene_level_TEProf2_truth",
    truth_label = "teprof2"
  )
})) |>
  mutate(depth = factor(depth, levels = depths)) |>
  arrange(depth)

saveRDS(teprof2_metrics_teinit, file.path(benchmark_dir, "TEProf2_official_by_depth_gene_level_metrics_TEProf2_truth.rds"))

# ------------------------------------------------------------------------------
# 7. Print Output
# ------------------------------------------------------------------------------
message("--- TEProf2 Evaluation Completed ---")
message("Metrics saved to: ", benchmark_dir)

