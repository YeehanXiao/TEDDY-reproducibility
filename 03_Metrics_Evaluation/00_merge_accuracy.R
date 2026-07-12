#!/usr/bin/env Rscript

# ============================================================
# Script: 00_merge_accuracy.R
# Cross-tool benchmark accuracy summary
#
# Input:
#   Precomputed gene-level benchmark metric RDS files for TEDDY,
#   Arriba, ChimeraTE, FREDY, LIONS, and TEProf2.
#
# Output:
#   1. official_benchmark_accuracy_summary.tsv
#   2. official_benchmark_accuracy_summary.rds
#
# Note:
#   Different tools were evaluated using task-aligned truth sets
#   corresponding to their intended output scope.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

# ------------------------------------------------------------
# 1. Parameters & Paths
# ------------------------------------------------------------
results_dir <- "./results"
summary_dir <- "./results/04_Summary"
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

depths <- c("5x", "10x", "25x", "50x", "100x")

tool_levels <- c(
  "TEDDY",
  "TEProf2",
  "FREDY",
  "LIONS",
  "Arriba",
  "ChimeraTE"
)

# ------------------------------------------------------------
# 2. Load precomputed metric tables
# ------------------------------------------------------------
# Ensure these paths perfectly match the outputs of the 03_Metrics_Evaluation scripts
metric_list <- list(
  TEDDY = readRDS(file.path(results_dir, "TEDDY_official_by_depth/TEDDY_official_by_depth_gene_level_metrics.rds")),
  
  Arriba = readRDS(file.path(results_dir, "arriba_merge/arriba_metrics_gene.rds")),
  
  ChimeraTE = readRDS(file.path(results_dir, "ChimeraTE_official_by_depth/evaluation/ChimeraTE_official_by_depth_gene_level_metrics.rds")),
  
  FREDY = readRDS(file.path(results_dir, "FREDY_official_by_depth/FREDY_official_by_depth_gene_level_metrics_FREDY50_truth.rds")),
  
  # Assuming LIONS metric extraction was done similarly to others
  LIONS = readRDS(file.path(results_dir, "LIONS_official_by_depth/LIONS_official_by_depth_gene_level_metrics_TEinitiated_truth.rds")),
  
  TEProf2 = readRDS(file.path(results_dir, "TEProf2_benchmark/TEProf2_official_by_depth_gene_level_metrics_TEProf2_truth.rds"))
)

# ------------------------------------------------------------
# 3. Define tool-specific benchmark scope and truth set
# ------------------------------------------------------------
tool_scope <- tibble(
  tool = c("TEDDY", "Arriba", "ChimeraTE", "FREDY", "LIONS", "TEProf2"),
  
  benchmark_scope = c(
    "TE-chimeric transcript reconstruction",
    "Fusion/breakpoint caller aggregated to gene level",
    "TE-gene chimeric event detection",
    "FREDY-compatible TE-overlap transcript detection",
    "Exon-repeat interaction / TE-initiated-compatible",
    "TE-initiated transcript annotation"
  ),
  
  truth_set = c(
    "Full expressed TE-chimeric gene",
    "Full expressed TE-chimeric gene",
    "Full expressed TE-chimeric gene",
    "FREDY-compatible >=50% TE-overlap gene",
    "Expressed TE-initiated gene",
    "Expressed TE-initiated gene"
  )
)

# ------------------------------------------------------------
# 4. Helper function to standardize metric columns
# ------------------------------------------------------------
extract_metric <- function(x, tool_name) {
  x <- as_tibble(x)
  
  expressed_col <- grep("^expressed_truth_gene", colnames(x), value = TRUE)
  
  if (length(expressed_col) == 0) {
    stop("No expressed_truth_gene column found for: ", tool_name)
  }
  
  pred_col <- dplyr::case_when(
    "pred_gene" %in% colnames(x) ~ "pred_gene",
    "chi_GTF_gene" %in% colnames(x) ~ "chi_GTF_gene",
    TRUE ~ NA_character_
  )
  
  if (is.na(pred_col)) {
    stop("No prediction gene column found for: ", tool_name)
  }
  
  x |>
    mutate(
      tool = tool_name,
      depth = as.character(depth),
      pred_gene_unified = .data[[pred_col]],
      expressed_truth_gene_unified = .data[[expressed_col[1]]]
    ) |>
    transmute(
      tool,
      depth,
      pred_gene = pred_gene_unified,
      expressed_truth_gene = expressed_truth_gene_unified,
      TP = .data[["TP"]],
      FP = .data[["FP"]],
      FN = .data[["FN"]],
      precision = .data[["precision"]],
      recall = .data[["recall"]],
      F1 = .data[["F1"]]
    )
}

# ------------------------------------------------------------
# 5. Build final accuracy summary table
# ------------------------------------------------------------
message("Aggregating metrics across all tools...")

accuracy_summary <- bind_rows(
  lapply(names(metric_list), function(nm) {
    extract_metric(metric_list[[nm]], nm)
  })
) |>
  left_join(tool_scope, by = "tool") |>
  mutate(
    depth = factor(depth, levels = depths),
    tool = factor(tool, levels = tool_levels)
  ) |>
  select(
    tool,
    benchmark_scope,
    truth_set,
    depth,
    pred_gene,
    expressed_truth_gene,
    TP,
    FP,
    FN,
    precision,
    recall,
    F1
  ) |>
  arrange(tool, depth)

# Export TSV and RDS
write_tsv(accuracy_summary, file.path(summary_dir, "official_benchmark_accuracy_summary.tsv"))
saveRDS(accuracy_summary, file.path(summary_dir, "official_benchmark_accuracy_summary.rds"))
