# ==============================================================================
# Script: 04_arriba_metrics.R
# Purpose: Aggregate and evaluate Arriba fusion calls at Gene & Region levels.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Paths and Constants (Desensitized)
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(rtracklayer)
  library(GenomicRanges)
})

truth_dir <- "./results"
arriba_dir <- "./results/arriba_merge"

truth_gtf_file <- file.path(truth_dir, "official_simulated_reference_90pct.gtf")
tx_truth_file <- file.path(truth_dir, "official_simulated_transcript_truth_status.tsv")
gene_truth_file <- file.path(truth_dir, "official_simulated_gene_truth_status.tsv")
iso_truth_file <- file.path(truth_dir, "official_simulated_1000.isoforms.results")

depths <- c("5x", "10x", "25x", "50x", "100x")

# ------------------------------------------------------------------------------
# 1. Read Truth Tables
# ------------------------------------------------------------------------------
tx_truth <- read.delim(tx_truth_file)
gene_truth <- read.delim(gene_truth_file)
iso_truth <- read.delim(iso_truth_file, check.names = FALSE)

tx_truth_expr <- tx_truth |>
  left_join(iso_truth |> select(transcript_id, TPM, expected_count), by = "transcript_id") |>
  mutate(TPM = coalesce(TPM, 0), expected_count = coalesce(expected_count, 0), is_expressed_truth = TPM > 1)

truth_gene_expr <- tx_truth_expr |> filter(is_TE_chimeric_tx, is_expressed_truth) |> distinct(gene_name)
truth_gene_raw <- gene_truth |> filter(is_TE_chimeric_gene) |> distinct(gene_name)

truth_gene_expr_detail <- tx_truth_expr |>
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

truth_summary <- data.frame(
  truth_TE_gene_raw = nrow(truth_gene_raw),
  truth_TE_gene_expressed = nrow(truth_gene_expr),
  truth_TE_tx_raw = sum(tx_truth$is_TE_chimeric_tx),
  truth_TE_tx_expressed = sum(tx_truth_expr$is_TE_chimeric_tx & tx_truth_expr$is_expressed_truth)
)

# ------------------------------------------------------------------------------
# 2. Read Arriba Outputs
# ------------------------------------------------------------------------------
arriba_files <- list.files(arriba_dir, pattern = "merge_.*_arriba_fusions(\\.discarded)?\\.tsv$", full.names = TRUE)
stopifnot(length(arriba_files) > 0)

parse_depth <- function(x) {
  sub("^merge_", "", sub("_arriba_fusions(\\.discarded)?\\.tsv$", "", basename(x)))
}

read_one_arriba <- function(f) {
  x <- fread(f, sep = "\t", header = TRUE, data.table = TRUE)
  if (nrow(x) == 0L) return(data.table())
  
  setnames(x, sub("^#", "", names(x)))
  required_cols <- c("gene1", "gene2", "breakpoint1", "breakpoint2")
  if (!all(required_cols %in% names(x))) stop("Missing required Arriba columns in: ", f)
  
  if (!"filters" %in% names(x)) x[, filters := ""]
  if (!"confidence" %in% names(x)) x[, confidence := NA_character_]
  
  split_cols <- intersect(c("split_reads1", "split_reads2"), names(x))
  mate_cols <- intersect(c("discordant_mates"), names(x))
  
  x[, `:=`(depth = parse_depth(f), source_file = basename(f), source_type = ifelse(grepl("\\.discarded\\.tsv$", basename(f)), "discarded", "final"))]
  
  x[, split_support := 0]
  if (length(split_cols) > 0L) x[, split_support := rowSums(as.data.frame(lapply(.SD, as.numeric)), na.rm = TRUE), .SDcols = split_cols]
  
  x[, discordant_support := 0]
  if (length(mate_cols) > 0L) x[, discordant_support := rowSums(as.data.frame(lapply(.SD, as.numeric)), na.rm = TRUE), .SDcols = mate_cols]
  
  x[, support := split_support + discordant_support]
  x[, event_id := paste(depth, gene1, gene2, breakpoint1, breakpoint2, sep = "|")]
  x
}

arriba_all <- rbindlist(lapply(arriba_files, read_one_arriba), fill = TRUE)
arriba_all[, filters_clean := tolower(ifelse(is.na(filters), "", filters))]
arriba_all[, `:=`(
  keep_relaxed_main = support >= 1,
  keep_relaxed_no_readthrough = support >= 1 & !grepl("read_through", filters_clean),
  keep_strict_main = source_type == "final" & support >= 2,
  keep_strict_no_readthrough = source_type == "final" & support >= 2 & !grepl("read_through", filters_clean)
)]

arriba_candidate_summary <- arriba_all[, .N, by = .(depth, source_type)][order(factor(depth, levels = depths), source_type)]

# ------------------------------------------------------------------------------
# 3. Gene-Level Evaluation
# ------------------------------------------------------------------------------
all_sim_gene <- gene_truth |> distinct(gene_name) |> filter(!is.na(gene_name), gene_name != "")

make_arriba_pred_gene_by_depth <- function(keep_col) {
  as_tibble(arriba_all) |>
    filter(.data[[keep_col]]) |>
    select(depth, event_id, gene1, gene2, support, filters, confidence, source_type) |>
    pivot_longer(cols = c(gene1, gene2), names_to = "partner", values_to = "gene_name") |>
    mutate(gene_name = as.character(gene_name)) |>
    separate_rows(gene_name, sep = ",") |>
    mutate(gene_name = trimws(gene_name), gene_name = sub("\\(.*\\)$", "", gene_name)) |>
    filter(!is.na(gene_name), gene_name != "") |>
    inner_join(all_sim_gene, by = "gene_name") |>
    distinct(depth, gene_name)
}

eval_arriba_gene_one_depth <- function(pred_gene_by_depth, depth_i, setting_label) {
  pred_gene <- pred_gene_by_depth |> filter(depth == depth_i) |> distinct(gene_name)
  
  tp_gene <- inner_join(pred_gene, truth_gene_expr, by = "gene_name")
  fp_gene <- anti_join(pred_gene, truth_gene_expr, by = "gene_name")
  fn_gene <- anti_join(truth_gene_expr, pred_gene, by = "gene_name")
  
  fp_gene_detail <- fp_gene |> left_join(truth_gene_raw |> mutate(in_raw_TE_truth = TRUE), by = "gene_name") |> mutate(in_raw_TE_truth = coalesce(in_raw_TE_truth, FALSE))
  fn_gene_detail <- fn_gene |> left_join(truth_gene_expr_detail, by = "gene_name")
  
  precision <- nrow(tp_gene) / (nrow(tp_gene) + nrow(fp_gene) + 1e-9)
  recall <- nrow(tp_gene) / (nrow(tp_gene) + nrow(fn_gene) + 1e-9)
  F1 <- 2 * nrow(tp_gene) / (2 * nrow(tp_gene) + nrow(fp_gene) + nrow(fn_gene) + 1e-9)
  
  data.frame(
    setting = setting_label, depth = depth_i, pred_gene = nrow(pred_gene),
    expressed_truth_gene_Arriba = nrow(truth_gene_expr), raw_truth_TE_gene_Arriba = nrow(truth_gene_raw),
    TP = nrow(tp_gene), FP = nrow(fp_gene), FN = nrow(fn_gene), precision = precision, recall = recall, F1 = F1,
    FP_in_raw_TE_truth = sum(fp_gene_detail$in_raw_TE_truth), FP_not_in_raw_TE_truth = sum(!fp_gene_detail$in_raw_TE_truth),
    FN_median_max_TPM = ifelse(nrow(fn_gene_detail) > 0, median(fn_gene_detail$max_TPM, na.rm = TRUE), NA_real_),
    FN_median_max_expected_count = ifelse(nrow(fn_gene_detail) > 0, median(fn_gene_detail$max_expected_count, na.rm = TRUE), NA_real_),
    FN_max_TPM = ifelse(nrow(fn_gene_detail) > 0, max(fn_gene_detail$max_TPM, na.rm = TRUE), NA_real_),
    FN_max_expected_count = ifelse(nrow(fn_gene_detail) > 0, max(fn_gene_detail$max_expected_count, na.rm = TRUE), NA_real_)
  )
}

arriba_pred_gene_main <- make_arriba_pred_gene_by_depth("keep_relaxed_main")
arriba_pred_gene_no_readthrough <- make_arriba_pred_gene_by_depth("keep_relaxed_no_readthrough")

arriba_pred_gene_summary <- bind_rows(
  arriba_pred_gene_main |> count(depth, name = "pred_gene") |> mutate(setting = "support >= 1, including read-through"),
  arriba_pred_gene_no_readthrough |> count(depth, name = "pred_gene") |> mutate(setting = "support >= 1, excluding read-through")
) |> mutate(depth_factor = factor(depth, levels = depths)) |> arrange(setting, depth_factor) |> select(-depth_factor)

arriba_metrics_sensitivity <- bind_rows(
  bind_rows(lapply(depths, \(depth_i) eval_arriba_gene_one_depth(arriba_pred_gene_main, depth_i, "support >= 1, including read-through"))),
  bind_rows(lapply(depths, \(depth_i) eval_arriba_gene_one_depth(arriba_pred_gene_no_readthrough, depth_i, "support >= 1, excluding read-through")))
) |> mutate(depth_factor = factor(depth, levels = depths), setting = factor(setting, levels = c("support >= 1, including read-through", "support >= 1, excluding read-through"))) |> arrange(setting, depth_factor) |> select(-depth_factor)

arriba_metrics_gene <- arriba_metrics_sensitivity |> filter(setting == "support >= 1, including read-through") |> mutate(setting = as.character(setting))
saveRDS(arriba_metrics_gene, file = file.path(arriba_dir, "arriba_metrics_gene.rds"))

# ------------------------------------------------------------------------------
# 4. Region-Level Sanity Overlaps
# ------------------------------------------------------------------------------
truth_gtf <- import(truth_gtf_file)
truth_tx <- truth_gtf[truth_gtf$type == "transcript"]
truth_tx_keep <- tx_truth_expr |> filter(is_TE_chimeric_tx, is_expressed_truth) |> pull(transcript_id) |> unique()
truth_tx_TE_expr <- truth_tx[as.character(truth_tx$transcript_id) %in% truth_tx_keep]
truth_txids_TE_expr <- as.character(truth_tx_TE_expr$transcript_id)

parse_bp_gr <- function(bp) {
  bp <- as.character(bp); chr <- sub(":.*$", "", bp); pos <- suppressWarnings(as.integer(sub("^.*:", "", bp)))
  ok <- !is.na(chr) & chr != "" & !is.na(pos)
  GRanges(seqnames = chr[ok], ranges = IRanges(start = pos[ok], end = pos[Pos[ok]]), strand = "*", idx = which(ok))
}

arriba_all[, idx := .I]
grA <- parse_bp_gr(arriba_all$breakpoint1); grB <- parse_bp_gr(arriba_all$breakpoint2)
ovA <- findOverlaps(grA, truth_tx_TE_expr, ignore.strand = TRUE); ovB <- findOverlaps(grB, truth_tx_TE_expr, ignore.strand = TRUE)

arriba_all[, `:=`(txA_TE_expr = rep(list(character()), .N), txB_TE_expr = rep(list(character()), .N))]

if (length(ovA) > 0L) {
  txA_dt <- data.table(idx = grA$idx[queryHits(ovA)], transcript_id = truth_txids_TE_expr[subjectHits(ovA)])[, .(txA_new = list(unique(transcript_id))), by = idx]
  arriba_all[txA_dt, txA_TE_expr := i.txA_new, on = "idx"]
}
if (length(ovB) > 0L) {
  txB_dt := data.table(idx = grB$idx[queryHits(ovB)], transcript_id = truth_txids_TE_expr[subjectHits(ovB)])[, .(txB_new = list(unique(transcript_id))), by = idx]
  arriba_all[txB_dt, txB_TE_expr := i.txB_new, on = "idx"]
}

arriba_all[, truth_TE_expr_tx_union := Map(union, txA_TE_expr, txB_TE_expr)]
arriba_all[, is_candidate_region_hit := lengths(truth_TE_expr_tx_union) > 0L]

make_arriba_region_metric <- function(keep_col, setting_label) {
  out <- arriba_all[get(keep_col) == TRUE, .(pred_event = uniqueN(event_id), TP_region_event = uniqueN(event_id[is_candidate_region_hit]), hit_truth_TE_expr_tx = uniqueN(na.omit(unlist(truth_TE_expr_tx_union)))), by = depth]
  out <- merge(data.table(depth = depths), out, by = "depth", all.x = TRUE)
  cols_to_fill <- c("pred_event", "TP_region_event", "hit_truth_TE_expr_tx")
  out[, (cols_to_fill) := lapply(.SD, \(x) fifelse(is.na(x), 0L, x)), .SDcols = cols_to_fill]
  out[, `:=`(setting = setting_label, truth_TE_expr_tx = length(truth_tx_keep), FP_region_event = pmax(pred_event - TP_region_event, 0L), FN_TE_expr_tx = pmax(length(truth_tx_keep) - hit_truth_TE_expr_tx, 0L), precision_region = fifelse(pred_event == 0L, NA_real_, TP_region_event / pred_event), recall_region_tx = hit_truth_TE_expr_tx / length(truth_tx_keep))]
  out[, F1_region := fifelse(is.na(precision_region) | is.na(recall_region_tx) | precision_region + recall_region_tx == 0, NA_real_, 2 * precision_region * recall_region_tx / (precision_region + recall_region_tx))]
  out[]
}

arriba_candidate_region_metric <- rbindlist(list(make_arriba_region_metric("keep_relaxed_main", "support >= 1, including read-through"), make_arriba_region_metric("keep_relaxed_no_readthrough", "support >= 1, excluding read-through")), fill = TRUE)
arriba_candidate_region_metric[, depth_factor := factor(depth, levels = depths)]; setorder(arriba_candidate_region_metric, setting, depth_factor); arriba_candidate_region_metric[, depth_factor := NULL]

# ------------------------------------------------------------------------------
# 5. Export Reviewer-Facing Outputs
# ------------------------------------------------------------------------------
saveRDS(
  list(
    truth_summary = truth_summary, arriba_metrics_gene = arriba_metrics_gene, arriba_metrics_sensitivity = arriba_metrics_sensitivity,
    arriba_candidate_summary = arriba_candidate_summary, arriba_candidate_region_metric = arriba_candidate_region_metric,
    arriba_pred_gene_by_depth = arriba_pred_gene_main, arriba_pred_gene_by_depth_no_readthrough = arriba_pred_gene_no_readthrough,
    arriba_pred_gene_summary = arriba_pred_gene_summary
  ),
  file.path(arriba_dir, "arriba_metric_final_reviewer_facing.rds")
)

message("--- Arriba Evaluation Completed ---")
print(arriba_metrics_gene)